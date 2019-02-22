pro	rfl,file,w,f,xdr=xdr,par=par,header=header

;+
;	Reading Lars' package computed spectra: a header  +
;	data (a line with wavelengths and multiple lines with fluxes)
;
;	IN: file	- string		name of the input file
;
;	OUT: w		- fltarr		wavelength
;	     f		- fltarr		fluxes (Hlambda, as computed
;							by the package
;
;	KEYWORDS:	xdr		when on, the data in the file
;					are save into an xdr file (file+'.xdr')
;
;			par		strarray with the parameters
;
;			header	strarray with the header
;
;	NOTE: when available in the local directory, rfl reads the xdr version
;		of the file requested.
;
;	C. Allende Prieto, UT, July 2005
;		"            , MSSL/UCL, July 2009 -- added header keyword
;-

if n_params() lt 3 then begin
	print,'% RFL: use -- rfl,file,w,f[,xdr=xdr,par=par]'
	return
endif

;command='if (-r '+file+'.xdr ) echo 1'
;spawn,command,back
back=file_test(file+'.xdr')
if (back) then begin
	file=file+'.xdr'
	restore,file
endif else begin
	h='hola'
	openr,10,file
	i=0
	while strpos(h,'* Data') lt 0 do begin
		readf,10,h
		if i eq 0 then header=h else header=[header,h]
		i=i+1
	endwhile
	readf,10,h
	w=strsplit(h,/extract)
	w=double(w[5:n_elements(w)-1])
	readf,10,h
	par=strmid(h,0,48)
	h=strmid(h,48,strlen(h))
	f=strsplit(h,/extract)
	f=double(f)
	d=[[w],[f]]
	while not eof(10) do begin
		readf,10,h
		sss=strmid(h,0,48)
		h=strmid(h,48,strlen(h))
		f=strsplit(h,/extract)
		f=double(f)
		d=[[d],[f]]		
		par=[par,sss]
	endwhile
	close,10	
endelse

if keyword_set(xdr) then begin
	save,d,par,header,file=file+'.xdr'
endif

ww=where(d[*,0] gt 0.0d0)
w=d[ww,0]
f=d[ww,1:n_elements(d[0,*])-1]

end
