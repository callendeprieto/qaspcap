pro	rd_apstar,filename,x,y,ey,bit=bit,hd=hd,sky=sky,esky=esky,$
	telluric=telluric,etelluric=etelluric,wcoef=wcoef,lsfcoef=lsfcoef,$
	vrad=vrad,vhelio=vhelio,verr=verr,snr=snr,ra=ra,dec=dec
	
;+
;	Reads an APOGEE apStar FITS file
;
;	IN:	filename	-string		Name of the file to read
;		
;	OUT: x		- dblarr		wavelength array (npix) [A]
;		 y		- fltarr		flux array (npix x nfiber) [1e-17 erg/cm2/s/A]
;		ey 		- fltarr		error array (npix x nfiber) [1e-17 erg/cm2/s/A];
;
;	KEYWORDS: bit	-intarr		Flags
;			  hd	-strarr		header
;			  sky	-fltarr		Sky flux (same units as y)
;			 esky   -fltarr		Error in the sky flux
;		 telluric   -fltarr		Telluric spectrum 
;		etelluric   -fltarr		Error in the telluric spectrum
;		    wcoef   -fltarr		Wavelength calibration coefficients
;		    lsfcoef   -fltarr		LSF coefficients
;			vrad	-float		Doppler shift (km/s)
;			vhelio	-float		Heliocentric radial velocity (km/s)
;			verr	-float		Estimated uncertainty in vhelio 
;			snr		-float		Est. S/N ratio
;			ra      -float		RA (deg)
;			dec     -float		DEC (deg)
;
; C. Allende Prieto, IAC, July 2011 (rd_apvisit)
; C. Allende Prieto, IAC, May 2012, adapted from rd_apvisit
;		"  , Radazul, August 2012, added ra/dec keywords
;-

if N_params() lt 3 then begin 
	print,'% rd_apstar: use -- rd_apstar,filename,x,y,ey[,bit=bit,hd=hd,sky=sky,esky=esky,telluric=telluric,etelluric=etelluric,wcoef=wcoef,vrad=vrad,vhelio=vhelio,verr=verr,snr=snr,ra=ra,dec=dec,lsfcoef=lsfcoef]'
	return
endif

y=mrdfits(filename,0,hd)
y=mrdfits(filename,1,hd1)
ey=mrdfits(filename,2,hd2)
if arg_present(bit) then bit=mrdfits(filename,3,hd3)
;x=mrdfits(filename,4,hd4)
if arg_present(sky) then sky=mrdfits(filename,4,hd4)
if arg_present(esky) then esky=mrdfits(filename,5,hd5)
if arg_present(telluric) then telluric=mrdfits(filename,6,hd6)
if arg_present(etelluric) then etelluric=mrdfits(filename,7,hd7)
if arg_present(lsfcoef) then lsfcoef=mrdfits(filename,8,hd8)

npix=sxpar(hd1,'NAXIS1')
wcoef=[sxpar(hd,'CRVAL1'),sxpar(hd,'CDELT1')]
x=10.d0^(dindgen(npix)*wcoef[1]+wcoef[0])

if arg_present(vrad) then vrad=sxpar(hd,'VRAD')
if arg_present(vhelio) then vhelio=sxpar(hd,'VHELIO')
if arg_present(verr) then verr=sxpar(hd,'VERR')
if arg_present(snr) then snr=sxpar(hd,'SNR')
if arg_present(ra) then ra=sxpar(hd,'RA')
if arg_present(dec) then dec=sxpar(hd,'DEC')

end
