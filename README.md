# WRF_boundary_coupling

Fortran code to couple global model passive tracer fields to WRF boundaries
Please read the general Readme file first:
psu-inversion/WRF_boundary_coupling/README.md



**Description of the global systems used as input fields for the WRF boundaries**

The code has been developed to couple model outputs from two global models:

# 1. GEOS-Chem global model
GEOS-Chem is a global 3-D chemical transport model (CTM) for atmospheric composition driven by meteorological input from the Goddard Earth Observing System (GEOS) of the NASA Global Modeling and Assimilation Office. It is applied by research groups around the world to a wide range of atmospheric composition problems. Scientific direction of the model is provided by the international GEOS-Chem Steering Committee and by User Working Groups. The model is managed by the GEOS-Chem Support Team, based at Harvard University and Dalhousie University with support from the US NASA Earth Science Division and the Canadian National and Engineering Research Council.

GEOS-Chem is a grass-roots community model owned by its users, and ownership implies some responsibilities as listed in our welcome page for new users. If you are interested in using GEOS-Chem, please contact the GEOS-Chem Support Team who will send you instructions for joining the user community and accessing the code.

To download the code sources and users' guides, visit the webpage: 
http://acmg.seas.harvard.edu/geos/

The version of the wrfchembc code for GEOS-Chem is in the directory "psu-inversion/WRF_boundary_coupling/wrfchem_GEOS_pkg/"
Main source code "main_bc_wrfchem_GEOS.f90"
Please refer to README_GEOSChem.md for more instructions.

# 2. CarbonTracker Inversion System (TM5 model)
CarbonTracker is a CO2 measurement and modeling system developed by NOAA to keep track of sources (emissions to the atmosphere) and sinks (removal from the atmosphere) of carbon dioxide around the world. CarbonTracker uses atmospheric CO2 observations from a host of collaborators and simulated atmospheric transport to estimate these surface fluxes of CO2. The current release of CarbonTracker, CT2016, provides global estimates of surface-atmosphere fluxes of CO2 from January 2000 through December 2015.

To represent the atmospheric transport, we use the Transport Model 5 (TM5). This is a community-supported model whose development is shared among many scientific groups with different areas of expertise. The model is developed and maintained jointly by the Institute for Marine and Atmospheric Research Utrecht (IMAU, The Netherlands), the Joint Research Centre (JRC, Italy), the Royal Netherlands Meteorological Institute (KNMI), the Netherlands Institude for Space Research (SRON), and the NOAA Earth System Research Laboratory (ESRL). 

For more information on the TM5 model, pleas visit:
http://tm5.sourceforge.net/

For more information on the CarbonTracker modeling system, please visit:
https://www.esrl.noaa.gov/gmd/ccgg/carbontracker/index.php

To access the original CarbonTracker model outputs, please visit:
https://www.esrl.noaa.gov/gmd/ccgg/carbontracker/molefractions.php

The version of the wrfchembc code for CarbonTracker is in the directory "psu-inversion/WRF_boundary_coupling/wrfchem_CT_pkg/"
Main source code "main_bc_wrfchem_CT.f90"
Please refer to README_GEOSChem.md for more instructions.

