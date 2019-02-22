pro	qwin,synthfile,csv,flt=flt

;+
;	Converts/adapts one or more window files in csv format for a particular
;	synthfile, creating FERRE-ready filter files (.flt).
;
;	INPUT: synthfile -string-	Name of the FERRE synthfile
;		csv	 -string- 	Name of the input csv file(s) 
;					(can work with wildcards, a string, 
;					or an array of strings)
;
;	KEYWORKDS: flt	 -string- 	Name of the output flt file
;					(if not given, set to the same
;					as the input csv file, suplemented
;					with the synthname and  ending in .flt, and
;					returned)
;
;	C. Allende Prieto, March 2016
;-

if N_Params() lt 2 then begin
	print,'% QWIN: -- usage -- qwin,synthfile,csv[,flt]'
	return
endif

if n_elements(csv) eq 1 then begin
	csv=file_search(csv)
	if max(where(csv ne '')) lt 0 then begin
		print,'% QWIN: error -- cannot find input files!'
		return
	endif
endif

if not keyword_set(flt) then flt=strarr(n_elements(csv))

;read synthfile wavelengths
read_synth,synthfile,/g,lambda=x
root=strmid(synthfile,2,strlen(synthfile)-5)

get_lun,lun
for i=0,n_elements(csv)-1 do begin

	;read csv file
	load,csv[i],d,/csv
	d=double(d)
	nel=n_elements(d[0,*])

	;define flt name
	flt[i]=strmid(csv[i],0,strlen(csv[i])-4)+'_'+root+'flt'

	openw,lun,flt[i]
	y=interpol([0.0,0.0,transpose(d[1,*]),0.0,0.0],$
		[0.0,d[0,0]-1d-6,transpose(d[0,*]),d[0,nel-1]+1d-6,1d10],x)
	printf,lun,transpose(y)
	close,lun


	if max(y) le 0.0 then $
	print,'% QWIN:  Warning -- no data different from zero in output file '+flt[i]

endfor
free_lun,lun

end


