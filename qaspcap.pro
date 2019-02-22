pro qaspcap,synthfile,indir,snr=snr,lsf2d=lsf2d,elem=elem,intype=intype,$
  condor=condor,slurm=slurm,ferrex=ferrex,contpar=contpar

;+
;	Quick ASPCAP 
;
;	IN: synthfile - string		Name of the synth database file
;		indir	  - string		Name of the directory with the input 
;						(apvisit/apstar/cdat/2column) files
;		
;	KEYWORDS: snr - float		When larger than 0, this value is used
;					to inject additional Gaussian noise in the
;					spectra, when they are prepared for ferre.
;
;		  lsf2d - integer	By default model fluxes in synthfile are
;	 				convolved with the instrumental kernel 
;					(the same for all fibers), but when lsf2d 
;					is switched on, custom LSF for each
;					fiber are read from apStarLSF files.
;
;		elem	- string	Filter files (.csv) for computing elemental
;					abundances. Wildcards work.
;
;		intype - string		Type of file (default is apStar FITS files)
;
;					apvisit -- apVisit FITS files from APOGEE 
;						(vacuum wavelengths)
;					apstar  -- apStar FITS files from APOGEE (vacuum)
;					cdat -- smoothed model files produced by 
;						iqconv (air)
;					2column -- two column ASCII files with 
;						spectra (air)
;					3column -- three column ASCII files with 
;						spectra (air)
;						(wavelength, flux, flux error)
;
;		condor	- 		prepares input files for condor
;
;		slurm	- 		prepares slurm script
;
;		ferrex  - 		FERRE executable (default is ~/ferre/src/a.out)
;
;		contpar -               Parameters for the continuum normalization --
;					when provided through this keyword they override
;					those specified in the synthfile
;
;	C. Allende Prieto, IAC, May-June 2012
;			 , Radazul, September 2012 -- added snr keyword
;			 , IAC, May 9, 2013 --  added SLURM script
;			 , IAC, Oct 17, 2014 -- added lsf2d
;			 , IAC, March 10, 2016 -- transform snr to float before using it
;			 , IAC, April 9, 2016 -- added keywords elem and intype
;			 , IAC, May 4, 2016 -- added 3column data type
;			 , La Palma, July 2016 -- added contpar and code taking into account mask
;-

if N_params() lt 2 then begin
	print,'% QASPCAP:  use -- qaspcap,synthfile,indir[,snr=snr,lsf2d=lsf2d,elem=elem,intype=intype,contpar=contpar] '
	return
endif

;environmental info
hostname=getenv('HOST')
pwd=getenv('PWD')
user=getenv('USER')

;ferre executable
if not keyword_set(ferrex) then ferrex='~/ferre/src/a.out'

;constants
clight=299792.458d0 ;km/s
sigma_to_fwhm=2.0*sqrt(-2.0*alog(0.5))
zlimit=0.01 ; ~ 3000 km/s --  accept 'stars' with |z|<zlimit

;extract wavelength info from library
read_synth,synthfile,/grid,hd=hd,npix=npix,ntot=ntot,ndim=ndim,ax=ax,$
        continuum=cont,wave=wave,logw=logw,lambda=lambda,nnpix=nnpix,vacuum=vacuum,resolution=resolution

;user-provided continuum parameters override those from te grid
if n_elements(contpar) gt 0 then cont=contpar

;read plate and zbest data files
if not keyword_set(intype) then intype='apstar'

case intype of
 'apvisit':   infiles=file_search(indir+'/apVisit-*',count=nfiles)
 'apstar':    infiles=file_search(indir+'/apStar-*',count=nfiles)
 'cdat':      infiles=file_search(indir+'/flux*cdat',count=nfiles)
 '2column':   infiles=file_search(indir+'/*',count=nfiles)
 '3column':   infiles=file_search(indir+'/*',count=nfiles)
else: intype='none'
endcase

if intype eq 'none' or nfiles lt 1 then begin
  if nfiles lt 1 then print,'% QASPCAP: -- error -- no valid files found' else $
    print,'% QASPCAP: -- error -- invalid intype option (must be apstar/apvisit/cdat/2column/3column)'
  return
endif

tmp=strsplit(indir,'/',/extract)
root=strcompress(tmp[0],/rem)
sroot=strcompress(strmid(synthfile,2,strlen(synthfile)-2)+'.'+root,/rem)


;aux files
nml=sroot+'.nml'
bas=sroot+'.bash'     ;bash script
con=sroot+'.submit'  ;condor script
slm=sroot+'.slurm'   ;SLURM script
inp=sroot+'.in.tar'  ;input data (frd/vrd/err)
out=sroot+'.out.tar' ;output data (opf?/mdl?)

;open output files
openw,1,sroot+'.vrd'
openw,2,sroot+'.frd'
openw,3,sroot+'.err'

winter=0
if keyword_set(lsf2d) then begin
  if intype ne 'apstar' and intype ne 'apvisit' then begin
    print,'% QASPCAP: Error -- lsf2d can only be used with apStar/apVisit files!'
    close,/all
    return
  endif
  openw,4,sroot+'.lsf'
  openw,5,sroot+'.wav'
  winter=2
endif

;loop over fibers in plate
print,'% QASPCAP: Processing fiber ...'

ic=0 ; counter tracking the actual numbre of admitted fibers
for i=0,nfiles-1 do begin
;for i=0,3-1 do begin

	print,'Reading spectrum #',i,' --',infiles[i],'...'
	
	case intype of 
	'apvisit': begin
		rd_apvisit,infiles[i],x,y,ey,bit=bit,hd=hd,wcoef=wcoef,$
			;sky=sky,esky=esky,telluric=telluric,etelluric=etelluric,$
			vrad=vrad,verr=verr,ra=ra,dec=dec,vhelio=vhelio,lsfcoef=lsfcoef
	   end 
	'cdat': begin
		rfl,infiles[i],x,y
		;load,infiles[i],dd
		;x=transpose(float(dd[0,*]))
		;y=transpose(float(dd[1,*]))
		ey=y/1000.d0
		bit=replicate(0,n_elements(y))
		vrad=0.0
		verr=0.0
		vhelio=0.0
		ra=0.0
		dec=0.0
	   end
	'2column': begin
		load,infiles[i],dd
		x=transpose(float(dd[0,*]))
		;airtovac,x
		y=transpose(float(dd[1,*]))
		ey=y/1000.d0
		bit=replicate(0,n_elements(y))
		vrad=0.0
		verr=0.0
		vhelio=0.0
		ra=0.0
		dec=0.0
	   end	   
	'3column': begin
		load,infiles[i],dd
		x=transpose(float(dd[0,*]))
		;airtovac,x
		y=transpose(float(dd[1,*]))
		ey=transpose(float(dd[2,*]))
		bit=replicate(0,n_elements(y))
		vrad=0.0
		verr=0.0
		vhelio=0.0
		ra=0.0
		dec=0.0
	   end	   	   
	'apstar': begin	
		rd_apstar,infiles[i],x,y,ey,bit=bit,hd=hd,wcoef=wcoef,$,
		;sky=sky,esky=esky,telluric=telluric,etelluric=etelluric,$
		vrad=vrad,verr=verr,ra=ra,dec=dec,vhelio=vhelio,lsfcoef=lsfcoef

		;retain only combined spectrum in the apStar files
		y=transpose(y[*,0])
		ey=transpose(ey[*,0])
		bit=transpose(bit[*,0])

	   end
	else: begin
			print,'% qaspcap: ERROR -- apvisit must be 0, 1 or 2'
			return
		end
	endcase
	
	;convert to vacuum if needed
	if min(vacuum) lt 1 and (intype eq 'apvisit' or intype eq 'apstar') then vactoair,x
	if min(vacuum) gt 0 and (intype eq 'cdat') then airtovac,x
	
	
	;find last . in filename 
	lastdotpos=strpos(infiles[i],'.',/reverse_search)
	object=strmid(infiles[i],0,lastdotpos)
	;replace '/' with '_'
	slashposdot=strpos(infiles[i],'/')
	for j=0,n_elements(slashposdot)-1 do strput,object,'_',slashposdot[j]
		
	z=vrad/clight
	
	;correct radial velocity and resample
	x2=x*(1.d0-z)
	
	
	if keyword_set(lsf2d) then begin
	    nlsf=145
	    xlambda=interpol(dindgen(n_elements(x)),x,lambda)
	    xstep=median(xlambda-shift(xlambda,1))
	    xlsf=dindgen(nlsf)*xstep-(nlsf-1)/2*xstep
	    xx=dblarr(npix,nlsf)
	    for j=0l,npix-1 do xx[j,*]=xlambda[j]+xlsf
	    lsfarr=lsf_gh2d(xx,xlambda,lsfcoef,dlsfgh)
	    wneg=where(lsfarr lt 0.0000d0)
	    if max(wneg) gt -1.0 then begin
	      if max(abs(lsfarr[wneg])) gt 1e-5 then stop,'% QASPCAP: Negative LSF!'
	      lsfarr[wneg]=0.0000d0
	    endif
	    for j=0l,npix-1 do lsfarr[j,*]=lsfarr[j,*]/total(lsfarr[j,*])
	endif else nlsf=0
	
	;stop
	
	if winter eq 2 then begin
	  nnpix2=nnpix
	  for j=0,n_elements(nnpix)-1 do begin
	    if j eq 0 then pix1=0 else pix1=long(total(nnpix[0:j-1]))
	    pix2=long(total(nnpix[0:j]))-1
	    wx2=where(x2 gt lambda[pix1+nlsf] and x2 lt lambda[pix2-nlsf] and y gt 1e-11)
	    nnpix2[j]=n_elements(wx2) ; will need this if continuum normalizing
	    print,'segment '+string(j+1),':',min(x2[wx2]),max(x2[wx2]),lambda[pix1+nlsf],lambda[pix2-nlsf]
	    if j eq 0 then begin
	      fluxarr=y[wx2]
	      efluxarr=ey[wx2]	  	  
	      bitarr=bit[wx2]
	      xarr=x2[wx2]
	    endif else begin
	      fluxarr=[fluxarr,y[wx2]]
	      efluxarr=[efluxarr,ey[wx2]]
	      bitarr=[bitarr,bit[wx2]]
	      xarr=[xarr,x2[wx2]]
	    endelse
	  endfor
	  
	  y=fluxarr
	  ey=efluxarr
	  bit=bitarr
	  x2=xarr
	    
	  print,'global ->',min(x2),max(x2),min(lambda),max(lambda)
	  
	endif else begin
	  nnpix2=nnpix
          wx2=where(y gt 1e-11)
	  y=interpol(y[wx2],x2[wx2],lambda)
	  ey=interpol(ey[wx2],x2[wx2],lambda)
	  bit=interpol(bit[wx2],x2[wx2],lambda)
	endelse
	
	if keyword_set(snr) then begin
		snr=float(snr) ; users may feed an integer
		y=y*(1.+randomn(seed,n_elements(y))/snr)
		ey=sqrt(ey^2+y^2/snr^2)
	endif
	
	
	;normalize
	if n_elements(cont) gt 1 then begin
	;if 1 eq 0 then begin
	  for j=0,n_elements(nnpix2)-1 do begin   
		if j eq 0 then pix1=0 else pix1=long(total(nnpix2[0:j-1]))
		pix2=long(total(nnpix2[0:j]))-1
		fluxarr=y[pix1:pix2]
		efluxarr=ey[pix1:pix2]
		bitarr=bit[pix1:pix2]

		mask=replicate(1,n_elements(bitarr))
		wbad=where(bitarr+1e-3 gt median(bitarr+1e-3,300))
		if max(wbad) gt -1 then mask[wbad]=0

		scont=size(cont)
		if scont[0] eq 2 then begin
		  continuum,cont[0,j],cont[1,j],cont[2,j],cont[3,j],fluxarr,fluxarr_cont,$
			mask=mask,ivar=1.d0/efluxarr^2
		  ;continuum,cont[0,j],cont[1,j],cont[2,j],cont[3,j],fluxarr,fluxarr_cont

		endif else begin
		  continuum,cont[0],cont[1],cont[2],cont[3],fluxarr,fluxarr_cont,$
			mask=mask,ivar=1.d0/efluxarr^2
		endelse

		;redo removing wild pixels (flux/cont >= 1.1 or flux/cont <=0) when cont[0,j] ge 0
		;if (cont[0,j] ge 0.0) then begin
		;  mask=intarr(pix2-pix1+1)
		;  wmask=where(fluxarr/fluxarr_cont lt 1.1 and $
		;	fluxarr/fluxarr_cont gt 0.00 $
		;	and efluxarr/fluxarr lt 1.00)
		;  if max(wmask) gt -1 then mask[wmask]=1
		;  continuum,cont[0,j],cont[1,j],cont[2,j],cont[3,j],fluxarr,fluxarr_cont,mask=mask
		;endif
		
		if j eq 0 then begin
		  fluxarr2=fluxarr/fluxarr_cont
		  efluxarr2=efluxarr/fluxarr_cont
		endif else begin
		  fluxarr2=[fluxarr2,fluxarr/fluxarr_cont]
		  efluxarr2=[efluxarr2,efluxarr/fluxarr_cont]
	        endelse	  
	          
	  endfor
	  y=fluxarr2
	  ey=efluxarr2

	endif 
	
	;plot,lambda,y
	;oplot,lambda,ey
	;stop
	
	;loadct,12
	;plot,fluxarr3,yr=[0.0,1.2],xr=[4000,6000]
	;oplot,y,col=180
	;print,mean(fluxarr3),median(fluxarr3),stdev(fluxarr3)
	;stop


	;cleanup useless data
	wbad=where(y le 0.0 or finite(y) eq 0)
	if max(wbad) gt -1 then begin
		y[wbad]=0.000000
		if n_elements(wbad) gt n_elements(y)*0.5 then pr=0
	endif
	wbad=where(ey le 0.0 or finite(ey) eq 0)
	if max(wbad) gt -1 then begin
		y[wbad]=0.000000		
		ey[wbad]=median(y)*1e5	
		if n_elements(wbad) gt n_elements(y)*0.5 then pr=0				
	endif	

	
	;print out output
	printf,1,object,0.0,0.0,0.0,0.0,vhelio,vrad,verr,ra,dec,$
			format='(a40,4(1x,f8.4),5(1x,e15.4))'
	printf,2,y,format='(1000000e14.6)'
	printf,3,ey,format='(1000000e14.6)'

	
	if keyword_set(lsf2d) then begin
	  for j=0l,nlsf-1 do printf,4,lsfarr[*,j],format='(1000000e14.6)'
	endif
	
	if winter eq 2 then begin
	  printf,5,x2,format='(1000000e14.6)'
	endif
	
	ic=ic+1
			
endfor
print,''
close,/all

;write main.nml
print,'% QASPCAP: writing '+nml
openw,1,nml
printf,1," &LISTA"
printf,1," NDIM = "+string(ndim,format='(i)')
printf,1," NOV = "+string(ndim,format='(i)')
printf,1," INDV = "+strjoin(string(indgen(ndim)+1),' ')
printf,1," SYNTHFILE(1) = '"+pwd+'/'+synthfile+"'"
printf,1," PFILE = '"+sroot+".vrd'"
printf,1," FFILE = '"+sroot+".frd'"
printf,1," ERFILE = '"+sroot+".err'"
printf,1," OPFILE = '"+sroot+".spm'"
printf,1," OFFILE = '"+sroot+".mdl'"
printf,1," SFFILE = '"+sroot+".nrd'"
printf,1," ERRBAR=1"
if keyword_set(lsf2d) then begin
	printf,1," LSF=12"
	printf,1," NLSF="+string(nlsf,format='(i4)')
	printf,1," LSFFILE = '"+sroot+".lsf'"
endif
if winter eq 2 then begin
	printf,1," WINTER="+string(winter,format='(i4)')
	printf,1," WFILE = '"+sroot+".wav'"
endif
indini=replicate(string(2),ndim)
;printf,1," INDINI = 2 1 1 1 3 2"
if n_elements(indini) gt 4 then indini[0:2]='1'
if n_elements(indini) gt 5 then indini[3]='1'
if n_elements(indini) gt 6 then indini[4]='1'
printf,1," NRUNS=",long(product(float(indini)))
printf,1," INDINI="+strjoin(indini,' ')
printf,1," COVPRINT=1"
printf,1," INTER=3"
printf,1," NTHREADS=1"
printf,1," F_FORMAT= 0"
printf,1," F_ACCESS = 0"
printf,1," CONT=1"
printf,1," NCONT=4"
;printf,1," PCAPROJECT=0"
;printf,1," PCACHI=0"
printf,1," /"
close,1

;elements
if keyword_set(elem) then begin

	;pick up synth file labels and sort the elements
	label=ax.label
	for i=0,n_elements(label)-1 do begin
		tmp1=strsplit(label[i],' ',/extract)
		tmp2=replicate(i,n_elements(tmp1))
		if i eq 0 then begin
			label_single=tmp1
			label_index=tmp2
		endif else begin
			label_single=[label_single,tmp1]
			label_index=[label_index,tmp2]
		endelse
	endfor
	remchar,label_single,"'"
	remchar,label_single," "
	metals_index=where(strmid(label,0,8) eq "'METALS'" or strmid(label,0,7) eq "'[M/H]'" or $
	                   strmid(label,0,8) eq "'[Fe/H]'")

	if n_elements(metals_index) ne 1 then begin
		print,'% QASPCAP: error -- cannot identify the label for the METALS dimension'
		return
	endif


	;read main nml file
	nlines=file_lines(nml)
	inputnml=strarr(nlines)
	openr,1,nml
	readf,1,inputnml
	close,1

	;prepare filter files and build element subdirs
 	qwin,synthfile,elem,flt=flt
	symbol=strarr(n_elements(elem))
	for i=0,n_elements(elem)-1 do begin
		symbol_starts=max(strpos(elem[i],'/'))+1
		symbol_ends=min(strpos(elem[i],'.filt'))
		symbol[i]=strmid(elem[i],symbol_starts,symbol_ends-symbol_starts)
		file_mkdir,symbol[i]
		cd,symbol[i]

		openw,1,'input.nml'
		for j=0,n_elements(inputnml)-1 do begin

			if (strcompress(inputnml[j],/rem) eq '/') then begin
				printf,1," FILTERFILE='../"+flt[i]+"'"
			endif
			if (strpos(inputnml[j],'NOV') gt -1) then begin
				printf,1," NOV=1"
				continue
			endif
			if (strpos(inputnml[j],'NRUNS') gt -1) then continue
			if (strpos(inputnml[j],'INDINI') gt -1) then continue
			if (strpos(inputnml[j],' PFILE') gt -1 or $
				;strpos(inputnml[j],' SYNTHFILE') gt -1 or $
				strpos(inputnml[j],' FFILE') gt -1 or $
				strpos(inputnml[j],' ERFILE') gt -1) then begin
					tmp1=strmid(inputnml[j],0,strpos(inputnml[j],"'")+1)	
					tmp2='../'
					tmp3=strmid(inputnml[j],strpos(inputnml[j],"'")+1,$
						strpos(inputnml[j],".",/reverse_search)-strpos(inputnml[j],"'"))
					tmp4=strmid(inputnml[j],strpos(inputnml[j],".",/reverse_search)+1,3)
					if strpos(inputnml[j],' PFILE') gt -1 then tmp4='spm'		
					if strpos(inputnml[j],' OPFILE') gt -1 then begin
						tmp2=''
						tmp4='opf'
					endif

					printf,1, tmp1+tmp2+tmp3+tmp4+"'"
					continue
			endif
			if (strpos(inputnml[j],'INDV') gt -1) then begin
				
				;by default very metallicity to find abundances
				vary_index=metals_index

				;but see if there is a better dimension too use
				wid=where(symbol[i] eq label_single)
				if (n_elements(wid) eq 1 and max(wid) gt -1) then $
					vary_index=label_index[wid]

				printf,1," INDV="+string(vary_index+1,format='(i)')
				continue
			endif
			printf,1,inputnml[j]
		endfor
		close,1
		cd,'..'
	endfor
endif


if keyword_set(condor) then begin
  ;spawn,'tar cvf '+inp+' '+sroot+'.vrd '+ sroot+'.frd '+$
  ;	sroot+'.lsf '+sroot+'.wav '+sroot+'.err '+nml,back,backerr

  synthfile_info=file_info(synthfile)
  frdfile_info=file_info(sroot+'.frd')
  
  openw,1,con
  printf,1,'executable = '+bas
  printf,1,'universe = vanilla'
  printf,1,'getenv = True'
  ;printf,1,'Requirements = TotalDisk > 5000000 && Memory > 1970 && LocalExecute == TRUE && Arch == "X86_64" && Machine != "carro.iac.es" && Machine != "ella.iac.es" && Machine != "indico.iac.es" && Machine != "jedey.iac.es" && Machine != "trevina.iac.es" && Machine != "espia.iac.es" && Machine != "chaplin.iac.es"  && Machine != "viga.ll.iac.es" && Machine != "cameron.ll.iac.es" && Machine != "trevina.ll.iac.es" && Machine != "voto.ll.iac.es" && Machine != "capo.iac.es" && Machine != "gere.iac.es" && Machine != "china.iac.es" && Machine != "echar.iac.es" && Machine != "ibero.iac.es" && Machine != "igual.iac.es" && Machine != "tacto.iac.es" && Machine != "tango.iac.es" && Machine != "veta.iac.es" && Machine != "viena.iac.es" && Machine != "embate.iac.es" && Machine != "cerca.iac.es" && Machine != "tapiz.iac.es" && Machine != "epico.iac.es" && Machine != "elche.iac.es" && Machine != "inicio.iac.es"&& Machine != "beso.iac.es"'
  printf,1,'request_disk   = '+string(frdfile_info.size/1e9*10.,format='(f6.1)')+' GB'
  printf,1,'request_memory = '+string(synthfile_info.size/1e9*1.2,format='(f5.1)')+' GB'
  printf,1,'Rank = KFlops'
  ;printf,1,'transfer_input_files = '+inp
  ;printf,1,'transfer_output_files = '+out
  ;printf,1,'WhenToTransferOutput = ON_EXIT'
  printf,1,'next_job_start_delay = 60'
  ;printf,1,'concurrency_limits = callende:400'
  printf,1,'output = '+sroot+'_condor.out'
  printf,1,'error =  '+sroot+'_condor.err'
  printf,1,'log =    '+sroot+'_condor.log'
  printf,1,'notification =   Never'
  printf,1,'queue'
  close,1

endif

openw,1,bas
printf,1,'#!/bin/bash'
printf,1,'echo started at '
printf,1,'date'
;if keyword_set(condor) then printf,1,'tar xvf '+inp
for i=0,0 do begin
	printf,1,"cp "+string(nml)+" input.nml"  
	printf,1,'echo running ferre with input.nml ...'
	printf,1,'cat input.nml'
	printf,1,'time '+ferrex+' >& '+sroot+'.flog'+string(i+1,format='(i1)')
endfor
for i=0,n_elements(elem)-1 do begin
  printf,1,"cd "+symbol[i]
  printf,1,'echo running ferre with input.nml ...'
  printf,1,'cat input.nml'
  printf,1,'time '+ferrex+' >& '+sroot+'.flog'+string(i+1,format='(i1)')
  printf,1,"cd .."
endfor
;if keyword_set(condor) then printf,1,'tar cvf '+out+' '+sroot+'.opf? '+sroot+'.mdl? '+sroot+'.flog?'
printf,1,'echo done at '
printf,1,'date'
spawn,'chmod ugo+x '+bas
close,1

if keyword_set(slurm) then begin
  openw,1,slm
  printf,1,'#!/bin/bash'
  printf,1,'# @ job_name  ='+sroot
  printf,1,'# @ initialdir = '+pwd
  printf,1,'# @ output = '+slm+'_%j.out'
  printf,1,'# @ error = '+slm+'_%j.err'
  printf,1,'# @ cpus_per_task = 4'
  printf,1,'# @ tasks_per_node = 1'
  printf,1,'# @ total_tasks = 1'
  printf,1,'# @ wall_clock_limit = 48:00:00'
  printf,1,'srun time '+bas
  close,1
endif

end
