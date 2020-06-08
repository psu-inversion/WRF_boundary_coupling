ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS =

if GEOSCHEM
bin_PROGRAMS += wrfchembc_GEOS
endif

if CARBONTRACKER
bin_PROGRAMS += wrfchembc_CT
endif

AM_FCFLAGS = $(NETCDF4_FFLAGS)
AM_LDFLAGS = $(NETCDF4_LDFLAGS)
LDADD = $(NETCDF4_FLIBS) $(NETCDF4_LIBS)

wrfchembc_GEOS_SOURCES = module_GEOS_lib.f90 module_wrfchem_lib.f90 main_bc_wrfchem_GEOS.f90
wrfchembc_CT_SOURCES = module_CT_lib.f90 module_wrfchem_lib.f90 main_bc_wrfchem_CT.f90

dist_doc_DATA = README.md README_CT.md README_GEOSChem.md LICENSE

MOSTLYCLEANFILES = module_ct_lib.mod module_geos_lib.mod	\
    module_wrfchem_lib.mod
