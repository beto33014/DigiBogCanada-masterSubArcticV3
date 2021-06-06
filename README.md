# DigiBog SubArctic
Modification to DigiBog for subartic climates which considers growing season temperatures and water tables for plant litter production and decomposition. Additionally, this model version considers groundwater flux. Model is intended to opearte at a weekly time step. More information about model development can be found in Ramirez et al. (2022).

Fotran files: DigiBog_MAIN_lumped_nopft.f90, DigiBog_Hydro_lumped.f90

Visual Studio 2019 project file: digiBog3d.vfproj

DigiBog SubArctic input files for peatland in Quebec:
010_DigiBog_BB_IN_information.txt: model paramters
020_DigiBog_BB_IN_net_rain.txt: weekly net precpitation (cm/yr) = rainfall + snowmelt - evapotranspiration
040_DigiBog_BB_IN_column_status.txt:
050_DigiBog_BB_IN_baltitude.txt:
051_DigiBog_BB_IN_gs_temp.txt:
052_DigiBog_BB_IN_gsBegin.txt:
053_DigiBog_BB_IN_gsEnd.txt:
054_DigiBog_BB_IN_gwFlux.txt:

