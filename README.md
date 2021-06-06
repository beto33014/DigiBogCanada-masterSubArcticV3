# DigiBog SubArctic

Modification to DigiBog for subarctic climates which considers growing season temperatures and water tables for calculations of plant litter production and decomposition. This model version considers groundwater flux and is intended to operate at a weekly time step. Detailed information about model development can be found in Ramirez et al. (2022). Model input files for a transect in Canada are provided.

Description of files in repository:

### Fortran source code
  *DigiBog\_MAIN\_lumped\_nopft.f90:* main model, decomposition and production processes  
  *DigiBog\_Hydro\_lumped.f90:* hydraulic processes
<br>

### Development environment
  *digiBog3d.vfproj:* visual studio 2019 project
<br>
### Input files 
  *010\_DigiBog\_BB\_IN\_information.txt*: model parameters  
  *020\_DigiBog\_BB\_IN\_net\_rain.txt*: weekly net precpitation (cm/yr) = rainfall + snowmelt - evapotranspiration  
  *040\_DigiBog\_BB\_IN\_column\_status.txt*: designate active columns and boundary conditions  
  *050\_DigiBog\_BB\_IN\_baltitude.txt*: elevation of mineral soil layer (cm)   
  *051\_DigiBog\_BB\_IN\_gs\_temp.txt*: time series of mean growing season temperature (&deg;C)      
  *052\_DigiBog\_BB\_IN\_gsBegin.txt*: time series of week of the year when the growing season begins    
  *053\_DigiBog\_BB\_IN\_gsEnd.txt*: time series of week of the year when the growing season ends     
  *054\_DigiBog\_BB\_IN\_gwFlux.txt*:  
