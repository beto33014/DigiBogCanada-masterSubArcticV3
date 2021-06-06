# DigiBog SubArctic

Modification to DigiBog for subartic climates which considers growing season temperatures and water tables for plant litter production and decomposition. This model version considers groundwater flux and is intended to opearte at a weekly time step. Detailed information about model development can be found in Ramirez et al. (2022).

Description of files in repository:

### Fortran code
  *DigiBog\_MAIN\_lumped\_nopft.f90:*  
  *DigiBog\_Hydro\_lumped.f90:*
<br>

### Visual Studio 2019
  *digiBog3d.vfproj:*
<br>
### Input files 
  *010\_DigiBog\_BB\_IN\_information.txt*: model parameters  
  *020\_DigiBog\_BB\_IN\_net\_rain.txt*: weekly net precpitation (cm/yr) = rainfall + snowmelt - evapotranspiration  
  *040\_DigiBog\_BB\_IN\_column\_status.txt*:  
  *050\_DigiBog\_BB\_IN\_baltitude.txt*:  
  *051\_DigiBog\_BB\_IN\_gs\_temp.txt*:    
  *052\_DigiBog\_BB\_IN\_gsBegin.txt*:  
  *053\_DigiBog\_BB\_IN\_gsEnd.txt*:    
  *054\_DigiBog\_BB\_IN\_gwFlux.txt*:  
