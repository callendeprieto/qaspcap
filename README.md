# qaspcap

An IDL wrapper for using FERRE on APOGEE spectra

qaspcap is a simple IDL script that substitutes the entire ASPCAP pipeline for quick-and-dirty work. It's easy to setup (you just need IDL, FERRE, and an ASPCAP spectral library) and runs both on apStar and apVisit files.

**Requirements**

- As much RAM as taken by the APOGEE model grid used. The usual PCA-compressed grids are a few GB (about 6 GB for the one in the example below). If you don't have it, you can still do it using the grid on disk, but at a slower speed (see the ​FERRE manual).
- IDL (or GDL)

- A FORTRAN95 or later compiler (gfortan 4.6 or later does just fine).

**Getting started**

Download qaspcap and its IDL dependencies. Unpack and make them visible to IDL.

Download an appropriate model grid from the APOGEE repository. 

Download and compile FERRE (additional notes are included at the bottom of the page).

Convert the ASPCAP grid to binary running ascii2bin (provided with FERRE). To to this simply call the executable (ascii2bin) and enter the name of the ASCII grid file when inquired by the program (e.g. p_apsasGK_131216_lsfcombo5v6_w123.dat), then answer unf when given the choice between formatted or unformatted (fmt/unf). The conversion takes about 20 minutes (it has to read the ascii file in to convert).

Get the spectra you want to analyze into a subdirectory. The code will run on all fits files in a directory at once.

**Atmospheric parameters**

At this point you are ready to run qaspcap to produce all the input files for FERRE, and run FERRE itself. You create a folder 'data' in the directory where your ASPCAP model grid is stored, and copy that file to it. Then run qaspcap:

IDL> qaspcap,'p_apsasGK_131216_lsfcombo5v6_w123.dat','data' 

% Compiled module: QASPCAP.

% Compiled module: READ_SYNTH.

% READ_SYNTH: Warning -- This is an NPCA grid -- 

% READ_SYNTH: --> npix will be adopted from the first header

% READ_SYNTH: Grid id is 'p_apsasGK_131216_lsfcombo5v6_w123.dat'

% READ_SYNTH: --> NP= 5 9 5 9 7 11 11

% READ_SYNTH: --> NPIX = 900

% Compiled module: GETAA.

% READ_SYNTH: x not present -- grid read, fluxes are not interpolated

% QASPCAP: Processing fiber ...

Reading spectrum # 0 --data/apStar-r8-2M17492453-2524258.fits...

% Compiled module: RD_APSTAR.

% Compiled module: MRDFITS.


% Compiled module: FXPOSIT.

% Compiled module: MRD_HREAD.

% Compiled module: FXPAR.

% Compiled module: GETTOK.

% Compiled module: VALID_NUM.

MRDFITS: Null image, NAXIS=0

% Compiled module: FXMOVE.

% Compiled module: MRD_SKIP.

MRDFITS: Image array (8575,4) Type=Real*4

MRDFITS: Image array (8575,4) Type=Real*4

MRDFITS: Image array (8575,4) Type=Int*2

MRDFITS: Image array (26) Type=Real*8

% Compiled module: SXPAR.

% Compiled module: INTERPOL

% Compiled module: CONTINUUM.

% Compiled module: POLY_FIT.

% Compiled module: POLY.

% QASPCAP: writing apsasGK_131216_lsfcombo5v6_w123.dat.data.nml

Unless you have the FERRE executable at the default location expected by qaspcap ('~/ferre/src/a.out') you will want to use the 'ferrex' keyword to point to it. For apVisit files, you will need to specify intype='apvisit'.

Note that in qaspcap you have to point to the ASCII version of the stellar library (not the result of ascii2bin). The code will not read in the entire file, only the header of the library to get the proper wavelength grid when it normalizes and resamples the input spectrum in prep for FERRE. If you point to the binary file, you will get an error message relating to "NNPIX" being undefined, which means that the synthetic spectrum was not read in properly in read_synth subroutine.

On output, qaspcap produces a bash script to run ferre (apsasGK_131216_lsfcombo5v6_w123.dat.data.bash in this example). All of the output scripts will be named <grid name>[dot]<directory>[dot]<file extension>.

The default mode for FERRE uses the ASCII version of the model grid, but if you have successfully converted the grid to binary with 'ascii2bin' (see above), you will want to set F_FORMAT=1 in the input file (apsasGK_131216_lsfcombo5v6_w123.dat.data.nml). This has to be changed manually, but the name of the input file can stay the name of the ascii file.

Note that the input file will be copied to 'input.nml' before invoking FERRE; for details see the ​FERRE manual. Go ahead and try running the bash script.

The grid in the example has 7 dimensions: log(micro), [C/Fe], [N/Fe], [alpha/Fe], [Fe/H], logg and Teff. On an Intel Xeon it takes 30 seconds to load the binary version of the grid in memory (this is done once, regardless of the number of spectra you are processing), and then 3 min to search for the solution (4 searches, since NRUNS=4). You can check if the output matches what I get for this run by downloading my output files​. The output paramers go to the OPFILE (apsasGK_131216_lsfcombo5v6_w123.dat.data.spm):

id	log(micro)	[C/Fe]	[N/Fe]	[alpha/Fe]	[Fe/H]	logg	Teff
data_apStar-r8-2M17492453-2524258	3.7248E-01	1.5332E-01	1.3574E-01	3.0224E-01	-6.6045E-01	7.5856E-01	3.6510E+03	...
In IDL you could examine the fitting typing (using routines in the tar ball you installed earlier)

IDL> load,'apsasGK_131216_lsfcombo5v6_w123.dat.data.nrd',f

IDL> load,'apsasGK_131216_lsfcombo5v6_w123.dat.data.mdl',m

IDL> read_synth,'p_apsasGK_131216_lsfcombo5v6_w123.dat',/grid,lambda=x

IDL> plot,x,f,yr=[0,1.4],xr=[15900.,16200.] 

IDL> oplot,x,m,col=180


