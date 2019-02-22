pro 	gconv,x,y,fwhm,nx,ny,ppr=ppr,nsigma=nsigma,inter=inter,original=original

;+
;	Gaussian convolution
;
;	IN:	x	-fltarr		independent variable
;		y	-fltarr		array of data to smooth 
;		fwhm	-float		FWHM of the Gaussian kernel (pixels)
;	
;	OUT:	nx	-fltarr		new sampled x-axis
;		ny	-fltarr		smoothed data
;		
;	KEYWORDS: ppr	-float	pixels per resolution element 
;						(default  is ppr=3.)
;	
;		  nsigma -float		how far the convolution goes
;
;		  inter  		if set, the x-axis is interpolated
;					to force a constant step; if the step
;					already constant, 'inter' will have no
;					effect
;
;		  original		if set, the original sampling is kept
;					for the output arrays, unless inter is
;					also set and an interpolation is 
;					performed, in which case, the output
;					x-axis  will just have the same
;					number of points as the input x, but
;					resampled to constant step.
;	
;	EXAMPLES:  1) Smooth an input spectrum with a Gaussian of FWHM=2.0 AA
;		assuming the spectrum is given on a uniform (linearly sampled)
;		wavelength scale with a step of 0.1 AA.
;
;		IDL> fwhm=2.0/0.1
;		IDL> gconv,w,f,fwhm,w2,f2
;
;		   2) Smooth the same spectrum with a Gaussian of FWHM=50. km/s
;
;		IDL> step=(max(alog(w))-min(alog(w)))/n_elements(w)
;		IDL> fwhm=50./299792.458/step
;		IDL> gconv,alog(w),f,fwhm,w2,f2,/inter
;
;		Because of the keyword inter, the spectrum will be internally
;		interpolated to a step in ln (w) of 
;		step=(max(alog(w))-min(alog(w)))/n_elements(w)
;		Then the 50 km/s (or V/c = 0.000166782) is equivalent to
;		50./299792.458/step pixels.
;	
;	C. Allende Prieto, Univ. of Texas, March 2006
;	""               , MSSL/UCL, July 2009 -- changed to double, avoid psf_gaussian
;					 , MSSL/UCL, August 2009 -- modified to avoid rate=0
;                    , IAC, January 2010 -- changed indgen -> indgen(,/long)    
;					 , IAC, January 2011 -- bug fixed (checking linear or log x scale)	
;  					 , IAC, November 2012 -- changed ppr type from integer to float
;-

if (N_params() lt 4) then begin
	message,'use -- gconv,x,y,fwhm,nx,ny[,ppr=ppr,nsigma=nsigma',/cont
	message,'                    inter=inter,original=original]',/cont
	return
endif

;initialize
tol=1d-3
nx=0.0d0
ny=0.0d0
if keyword_set(inter) then begin
	nel=n_elements(x)
	xinput=x
	yinput=y
	minx=min(x)
	step=(max(x)-minx)/nel
	x=dindgen(nel)*step+minx
	y=interpol(yinput,xinput,x)
endif else begin
	;checking that the input x is evenly sampled
	nel=n_elements(x)				; nel input pixels
	step=abs(x-shift(x,1))
	med=median(step[1:nel-1])
	if max(abs(step[1:nel-1]-med)) gt tol then begin
		;not linearly spaced, maybe logarithmically
		lx=alog10(x)
		step=abs(lx-shift(lx,1))
		med=median(step[1:nel-1])
		if max(abs(step[1:nel-1]-med)) gt tol then begin
			message,'not log spaced, cannot handle it without interpolation',/cont
			return			
		endif
	endif
endelse

;set default params
if not keyword_set(ppr) then ppr=3.0     		;sample output function with ppr
						;points per FWHM		
if not keyword_set(nsigma) then nsigma=3.0	;how far the convolution goes
;Gaussian kernel				; units are pixels
sigma=fwhm/2.d0/sqrt(-2.d0*alog(0.5d0))		; equivalent sigma for kernel
grange=floor(nsigma*sigma)			; range for kernel (-range:range)
;psf=psf_gaussian(npix=2*grange+1,fwhm=fwhm,/normalize,ndimen=1)
psf=1.d0/sqrt(2.d0*!dpi)/sigma*exp(-((dindgen(2*grange+1)-grange)/sigma)^2/2.d0)
psf=psf/total(psf)

;select output x-axis
if keyword_set(original) then begin
;keep input sampling on output arrays
	sampling=indgen(nel,/long)
endif else begin
	;rate=floor(fwhm)/ppr            ; old version, whih ppr being integer
	rate=floor(fwhm/ppr)			; 1/rate pixels will be kept
	;make sure that rate is at least 1
	if rate lt 1 then rate=1
	sampling=uniq(indgen(nel,/long)/rate)-rate/2
endelse


;trim edges
wtrim=where(sampling ge grange[0] and sampling lt nel-grange[0])
if max(wtrim) le 0 then begin
	message,'The Gaussian is too wide for the input range',/cont
	return
endif
sampling=sampling[wtrim]
nx=x[sampling]					; new x-axis
nel2=n_elements(nx)				; nel2 output pixels
ny=dblarr(nel2)


;convolve
for i=0l,nel2-1 do begin
	ny[i]=total(psf*y[sampling[i]-grange:sampling[i]+grange])
	;print,'y=',y[sampling[i]-grange:sampling[i]+grange]
	;print,'psf=',psf
	;print,'ny[i]=',ny[i]
	;stop
endfor

;return the original values of x and y if resampled
if keyword_set(inter) then begin
	x=xinput
	y=yinput
endif

end
