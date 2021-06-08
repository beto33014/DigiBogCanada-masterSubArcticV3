PROGRAM DIGIBOG_3D
!New version as of 06/06/2021. Minor modifications tracked via the DigiBog 
!GitHub repository.

! -----------------------------------------------------------------------------
! Section 1.0 Program header
! -----------------------------------------------------------------------------

  !Description
  !A model to simulate water tables  and peat accumulation / wastage in
  !ombrotrophic bogs and other shallow aquifers using the Boussinesq equation.
  !The model uses centimetres and seconds.
  !Details of how the model works may be obtained from the code owners
  !(email addresses below).

  !Summary of code configuration
  !Net rainfall = annual, monthly or weekly
  !Temperature = annually
  !Layers = lumped
  !Recalcitrance = oxic and anoxic

  !Current code owners
  !Paul J. Morris, Andy J. Baird*, Lisa R. Belyea, Dylan M. Young, Pete Gill
  !*School of Geography
  !University of Leeds
  !Leeds
  !LS2 9JT
  !p.j.morris@reading.ac.uk
  !a.j.baird@leeds.ac.uk
  !l.belyea@qmul.ac.uk
  !geodmy@leeds.ac.uk
  !gy12pg@leeds.ac.uk


  !Modification history of code
  !Programmer           Date           Modifications
  !============================================================================
  !Andrew J. Baird      08/04/2005     Original 1.5-d code
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/02/2006     Conversion to 2.5-d and Fortran
  !                                    Standards implemented
  !----------------------------------------------------------------------------
  !Paul J. Morris       27/02/2006     Testing, minor corrections
  !----------------------------------------------------------------------------
  !Paul J. Morris       05/03/2006     Subroutines 'column_activation' and
  !                                    'steady_state_check' written. Former
  !                                    no longer exists (removed 18/07/2013).
  !----------------------------------------------------------------------------
  !Paul J. Morris       19/06/2006     Testing completed
  !----------------------------------------------------------------------------
  !Paul J. Morris       20/03/2007     Code cleaning
  !----------------------------------------------------------------------------
  !Paul J. Morris       24/04/2007     Further cleaning, including removal of
  !                                    replicated spatial step reference in
  !                                    'move_water' subroutine
  !----------------------------------------------------------------------------
  !Paul J. Morris       09/05/2007     Above-ground storage facilitated in
  !                                    2.5-d version
  !----------------------------------------------------------------------------
  !Paul J. Morris       02/07/2008     Final code cleaning, full annotation
  !----------------------------------------------------------------------------
  !Andy J. Baird        18/07/2013     Addition of variable net rainfall read
  !                                    in from file. Change to how boundary
  !                                    conditions specified; zero-flow
  !                                    Neumann condition also now allowed for,
  !                                    as are internal boundaries. Change to
  !                                    how output written to file (can now
  !                                    write to file multiple times before end
  !                                    of run).
  !----------------------------------------------------------------------------
  !Andy J. Baird        27/08/2013     Code cleaning and de-bugging of version
  !                                    from 18/07/2013
  !----------------------------------------------------------------------------
  !Dylan M. Young       20/05/14       Combine 1D and Hydro models for use in
  !                                    blanket bog simulation of management
  !                                    impacts. For testing, set fixed Dirichlet
  !                                    condition of column elevation * 0.5.
  !                                    Added memory allocation for x_flux and
  !                                    y_flux in wat_k_mean (Hydro_procs) for
  !                                    de-bugging purposes.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       11/06/14       Reorganisation of code as suggested by
  !                                    AJB and agreed with AJB and PJM.
  !                                    Removal of de-bugging code and resetting
  !                                    x_flux and y_flux to local arrays.
  !                                    Removal of steady-state check.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       18/06/14       Addition of variable Dirichlet condition
  !                                    water-table read in from the parameter
  !                                    file.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       13/07/14       Added timestep multiplier for calculated
  !                                    timestep for stability improvements
  !-----------------------------------------------------------------------------
  !Dylan M. Young       01/01/15       Added mineral layer for sloping model
  !                                    stability. Can be used for both raised
  !                                    and blanket bog models.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       17/03/15       Weather now read-in on a weekly basis
  !                                    (can be used for stable or variable
  !                                    weather). Mean annual temp. is used for
  !                                    new litter production.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       20/03/15       Recalcitrance added to anoxic
  !                                    decomposition.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       28/04/15       Added write-out for K profile and maximum
  !                                    timestep calc.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       06/12/16       Updated recalcitrance for anoxic
  !                                    decomp based on Morris et. al.,2015
  !                                    Rainfall varies weekly, temperature is
  !                                    set on an annual basis
  !-----------------------------------------------------------------------------
  !Pete Gill            04/2017        Added layer lumping routine to reduce
  !                                    model run time.
  !-----------------------------------------------------------------------------
  !Dylan M. Young       06/07/17       Updated recalcitrance to include oxic and
  !                                    anoxic
  !                                    decomp based on Morris et. al.,2015
  !-----------------------------------------------------------------------------
  !Dylan M. Young       13/05/20       New version including bug fixes in the 
  !                                    lumping routine for calculating the 
  !                                    layer remaining mass correctly and to 
  !                                    add the recalculation of transmissivity
  !                                    after layers have been lumped. This 
  !                                    version supersedes all previous versions
  !                                    of the ecohydrological model code.
  !-----------------------------------------------------------------------------
  !Jorge A. Ramirez     6/06/21        Modifications to DigiBog were made for 
  !                                    subarctic climates by modifying plant 
  !                                    litter production and peat decomposition.
  !                                    Groundwater flux as additional input. 
  !----------------------------------------------------------------------------- 
  !Code description.
  !Language: Fortran 95.
  !Software standards: Andrews, P. (1998). Unified Model Documentation Paper No.
  !3, Software Standards for the Unified Model: Fortran and Unix, Version 7.2,
  !Met Office, UK, 17 pp.
  !-----------------------------------------------------------------------------

  !Modules used:
  USE hydro_procedures_transmax
  IMPLICIT NONE
  1100 FORMAT (1X,I3,I3,I5,I5,4(20F20.10)) 

!-------------------------------------------------------------------------------
! Section 2.0 Definitions; Declarations
!-------------------------------------------------------------------------------

! Single-valued variable (and model parameter) declarations.

  INTEGER :: alloc_error,        &  !Memory allocation error flag
             t_extent,           &  !Number of simulated model years
             output_interval,    &  !Interval (years) for writing to output
                                    !files
             sub_year_counter,   &  !Counter for use in main model loop
             sub_year_counter_su,&  !Counter for use in main model loop
                                    !for summer water tables
             output_counter,     &  !Output counter for use in main time loop
             year_counter,       &  !Counter for total number of model years
             sub_annual_counter, &  !Counter for total weeks or months
             new_net_rain,       &  !Flag for reading in new net rain 
             x_extent,           &  !Model grid extent in x direction (number)
             y_extent,           &  !Model grid extent in y direction (number)
             z_extent,           &  !Spatial index counter
             x,                  &  !Spatial index counter
             y,                  &  !Spatial index counter
             z,                  &  !Spatial index counter
             net_rain_res,       &  !Net rainfall temporal resolution
             layer_exclusion,    &  !Don't lump these layers (numbered from top
                                    !of a column)
             growingSeasonBegin, &  !week growing season begins
             growingSeasonEnd,   &  !week growing season ends
             lump_index

  REAL(kind=q) :: oxic_decay_base,     &  !Proportion (yr^-1)
                  anoxic_decay_base,   &  !Proportion (yr^-1)
                  base_temp,           &  !deg. Celsius
                  Q10_oxic,            &  !Factor (dimensionless)
                  Q10_anoxic,          &  !Factor (dimensionless)
                  oxic_decay,          &  !Proportion (yr^-1)
                  anoxic_decay,        &  !Proportion (yr^-1)
                  density,             &  !Peat bulk density (g cm^-3)
                  porosity,            &  !Drainable porosity (cm^3 cm^-3 --
                                          !a proportion)
                  k_param_a,           &  !The a parameter in the k equation
                                          !(cm yr^-1)
                  k_param_b,           &  !The b parameter in the k equation
                                          !(dimensionless)
                  timestep,            &  !yr
                  t_min,               &  !min timestep
                  t_max,               &  !max timestep
                  t_min_years,         &  !No of years before timestep increases
                  t_max_years,         &  !Period over which timestep increases
                  t_rate,              &  !Rate of increase in timestep
                  transmax,            &  !Maximum transmissivity for timestep
                  t_step_multi,        &  !timestep multiplication factor
                  t_step_sum_netr,     &  !sum of sub annual timesteps
                  t_step_sum_su,       &  !sum of summer timesteps
                  mean_t_step,         &  !Average annual timesteps
                  annual_tsteps,       &  !Number of annual timesteps
                  net_rain_tstep,      &  !Number of weekly timesteps
                  rainfall,            &  !Net rainfall (cm yr^-1)
                  temperature,         &  !deg. C
                  gsTemperature,       &  !mean growing season temperature deg. C
                  spatial_step,        &  !Model horizontal increment
                                          !(x & y) (cm)
                  pond_depth,          &  !Depth of surface ponding (cm)
                  time_start,          &
                  time_finish,         &
                  time_elapsed,        &
                  weight_z,            &
                  weight_z_minus,      &
                  threshold,           & !Lump layers < this height (cm)
                  mineral,             & !mineral soil thickness (cm)
                  oxic_recal,          & !oxic recalcitrance exponent
                  anoxic_recal           !anoxic recalcitrance exponent

  !Climate intervals
  CHARACTER(LEN=1) :: weather_error          !error flag

  !Input and output files
  CHARACTER(LEN=40) :: data_file_name_010, & !IN Model run information file
                      data_file_name_020, &  !IN net rainfall values
                      data_file_name_040, &  !IN column activation status file
                      data_file_name_050, &  !IN base altitude input file
                      data_file_name_051, &  !IN mean growing season temperature per year
                      data_file_name_052, &  !IN week when growing season begins (1-52)
                      data_file_name_053, &  !IN week when growing season ends (1-52)
                      data_file_name_054, &  !In groundwater flux per column (cm/yr)
                      data_file_name_060, &  !OUT Layer mass remaining output
                      data_file_name_070, &  !OUT Column height output file
                      data_file_name_080, &  !OUT Water-table height output file
                      data_file_name_090, &  !OUT Transmissivity output file
                      data_file_name_100, &  !OUT Mass per area output file
                      data_file_name_110, &  !OUT Mean annual column water table
                      data_file_name_120, &  !OUT Layer wet proportion file
                      data_file_name_130, &  !OUT timestep
                      data_file_name_140, &  !OUT timestep sum
                      data_file_name_160, &  !OUT Summer water table
                      data_file_name_170, &  !OUT Run Time
                      data_file_name_180, &  !OUT Layers per column
                      data_file_name_190, &  !OUT Layer mass remaining, 
                                             !elevation, and thickness
                      data_file_name_200, & !OUT initial layer mass               
                      data_file_name_210, & !OUT final layer height
                      data_file_name_220, & !OUT layer age
                      data_file_name_230    !OUT annual column production

  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:) :: base_altitude,     & !Above datum
                                               ground_water_flux, & !cm/yr
                                               water_change,      &
                                               water_table,       & !Above base
                                               wk_mean,           & !Depth-av. K
                                               col_wt_depth,      & !cm
                                               col_wt_sum,        & !cm
                                               col_wt_sum_su,     & !cm
                                               col_wt_depth_av,   & !cm
                                               col_wt_depth_av_su,& !cm
                                               col_mass_per_area    !g cm^-2

  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: no_layers  !Number of layers (x,y)

  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: layers_in_layer 


  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:,:) :: wet_proportion, & !proportion
                                                 layer_age !, &
                                                 !layers_in_layer

  !layer_attributes stores layer thickness, k (cm yr^-1) and s (cm^3 cm^-3)
  !transmissivity stores layer elevation above base and transmissivity
  !layer_storage stores layer elevation (cm) above base and each layer's water
  !capacity as volume per unit area; i.e. expressed as a depth (cm)
  !layer_mass stores layer current (g cm^-2), initial (g cm^-2)
  !and remaining mass (proportion).
  REAL(KIND=q), ALLOCATABLE, DIMENSION(:,:,:,:) :: layer_attributes, &
                                                   transmissivity, &
                                                   layer_storage, &
                                                   layer_mass

  !activation_status indicates whether/how column participates in simulation
  CHARACTER(8), ALLOCATABLE, DIMENSION(:,:) :: activation_status

!-------------------------------------------------------------------------------
! Section 3.0 Data Input; Memory Management
!-------------------------------------------------------------------------------

  ! Name input and output data files.

  !Inputs
  data_file_name_010 = "010_DigiBog_BB_IN_information.txt"
  data_file_name_020 = "020_DigiBog_BB_IN_net_rain.txt"
  data_file_name_040 = "040_DigiBog_BB_IN_column_status.txt"
  data_file_name_050 = "050_DigiBog_BB_IN_baltitude.txt"
  data_file_name_051 = "051_DigiBog_BB_IN_gs_temp.txt"
  data_file_name_052 = "052_DigiBog_BB_IN_gsBegin.txt"
  data_file_name_053 = "053_DigiBog_BB_IN_gsEnd.txt"
  data_file_name_054 = "054_DigiBog_BB_IN_gwFlux.txt"
  
  !Outputs
  data_file_name_060 = "060_DigiBog_BB_OUT_layer_mass.txt"
  data_file_name_070 = "070_DigiBog_BB_OUT_column_height.txt"
  data_file_name_080 = "080_DigiBog_BB_OUT_wt_height.txt"
  data_file_name_090 = "090_DigiBog_BB_OUT_transmiss.txt"
  data_file_name_100 = "100_DigiBog_BB_OUT_mass_area.txt"
  data_file_name_110 = "110_DigiBog_BB_OUT_wt_depth.txt"
  data_file_name_120 = "120_DigiBog_BB_OUT_layer_wet_prop.txt"
  data_file_name_130 = "130_DigiBog_BB_OUT_col_t_step.txt"
  data_file_name_140 = "140_DigiBog_BB_OUT_t_step_sum.txt"
  data_file_name_160 = "160_DigiBog_BB_OUT_wt_depth_summer.txt"
  data_file_name_170 = "170_DigiBog_BB_OUT_run_time.txt"
  data_file_name_180 = "180_DigiBog_BB_OUT_layers_in_columns.txt"
  data_file_name_190 = "190_DigiBog_BB_OUT_rem_mass_file.txt"
  data_file_name_200 = "200_DigiBog_BB_OUT_init_mass_file.txt"
  data_file_name_210 = "210_DigiBog_BB_OUT_layer_ht_file.txt"
  data_file_name_220 = "220_DigiBog_BB_OUT_layer_age_file.txt"
  data_file_name_230 = "230_DigiBog_BB_OUT_column_production.txt"

! Open files for input and output of data.
  OPEN(UNIT=010, FILE=data_file_name_010, STATUS="OLD")
  OPEN(UNIT=020, FILE=data_file_name_020, STATUS="OLD")
  OPEN(UNIT=040, FILE=data_file_name_040, STATUS="OLD")
  OPEN(UNIT=050, FILE=data_file_name_050, STATUS="OLD")
  OPEN(UNIT=051, FILE=data_file_name_051, STATUS="OLD")
  OPEN(UNIT=052, FILE=data_file_name_052, STATUS="OLD")
  OPEN(UNIT=053, FILE=data_file_name_053, STATUS="OLD")
  OPEN(UNIT=054, FILE=data_file_name_054, STATUS="OLD")
  OPEN(UNIT=060, FILE=data_file_name_060, STATUS="REPLACE")
  OPEN(UNIT=070, FILE=data_file_name_070, STATUS="REPLACE")
  OPEN(UNIT=080, FILE=data_file_name_080, STATUS="REPLACE")
  OPEN(UNIT=090, FILE=data_file_name_090, STATUS="REPLACE")
  OPEN(UNIT=100, FILE=data_file_name_100, STATUS="REPLACE")
  OPEN(UNIT=110, FILE=data_file_name_110, STATUS="REPLACE")
  OPEN(UNIT=120, FILE=data_file_name_120, STATUS="REPLACE")
  OPEN(UNIT=130, FILE=data_file_name_130, STATUS="REPLACE")
  OPEN(UNIT=140, FILE=data_file_name_140, STATUS="REPLACE")
  OPEN(UNIT=160, FILE=data_file_name_160, STATUS="REPLACE")
  OPEN(UNIT=170, FILE=data_file_name_170, STATUS="REPLACE")
  OPEN(UNIT=180, FILE=data_file_name_180, STATUS="REPLACE")
  OPEN(UNIT=190, FILE=data_file_name_190, STATUS="REPLACE")
  OPEN(UNIT=200, FILE=data_file_name_200, STATUS="REPLACE")
  OPEN(UNIT=210, FILE=data_file_name_210, STATUS="REPLACE")
  OPEN(UNIT=220, FILE=data_file_name_220, STATUS="REPLACE")
  OPEN(UNIT=230, FILE=data_file_name_230, STATUS="REPLACE")

! Read data from information input data file.
  READ (010, *) oxic_decay_base
  READ (010, *) anoxic_decay_base
  READ (010, *) base_temp
  READ (010, *) Q10_oxic
  READ (010, *) Q10_anoxic
  READ (010, *) density
  READ (010, *) porosity
  READ (010, *) k_param_a
  READ (010, *) k_param_b
  READ (010, *) t_extent
  READ (010, *) annual_tsteps
  READ (010, *) output_interval
  READ (010, *) x_extent
  READ (010, *) y_extent
  READ (010, *) spatial_step
  READ (010, *) pond_depth
  READ (010, *) t_step_multi
  READ (010, *) layer_exclusion
  READ (010, *) threshold
  READ (010, *) oxic_recal
  READ (010, *) anoxic_recal
  READ (010, *) mineral
  READ (010, *) net_rain_res
  READ (010, *) t_min
  READ (010, *) t_max
  READ (010, *) t_min_years
  READ (010, *) t_max_years

  !Check if the climate resolution is correct
  WRITE (*,*) "Net rain resolution =",net_rain_res, "(1=year, 2=month, 3=week)"

  WRITE (*, '(A31)') "Is this  correct (Y/N)?"
  READ *, weather_error 
  SELECT CASE ( weather_error )
    CASE ("n")
      WRITE (*, *) "Incorrect interval for climate input"
        STOP
    CASE ("N")
      WRITE (*, *) "Incorrect interval for climate input"
        STOP
    CASE DEFAULT
      WRITE (*, *) "Climate resolution okay - simulation continues"
  END SELECT


  !Set the extent of layers
  z_extent = t_extent + 2 !  Ponding layer + added mineral + peat layers

  !Allocate memory to arrays.

  !Columns
  ALLOCATE(col_wt_depth(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth"
    STOP
  END IF

  ALLOCATE(col_wt_depth_av(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth_av"
    STOP
  END IF

  ALLOCATE(col_wt_depth_av_su(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_depth_av"
    STOP
  END IF

  ALLOCATE(col_wt_sum(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_sum"
    STOP
  END IF

  ALLOCATE(col_wt_sum_su(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_wt_sum_su"
    STOP
  END IF

  ALLOCATE(col_mass_per_area(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for col_mass_per_area"
    STOP
  END IF

  ALLOCATE(base_altitude(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for base_altitude"
    STOP
  END IF
  
  ALLOCATE(ground_water_flux(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for ground_water_flux"
    STOP
  END IF

  ALLOCATE(water_change(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_change"
    STOP
  END IF

  ALLOCATE(water_table(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for water_table"
    STOP
  END IF

  ALLOCATE(wk_mean(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for wk_mean"
    STOP
  END IF

  ALLOCATE(activation_status(x_extent, y_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for activation_status"
    STOP
  END IF

  ALLOCATE(wet_proportion(x_extent, y_extent, z_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for wet_proportion"
    STOP
  END IF

  ALLOCATE(layer_age(x_extent, y_extent, z_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for wet_proportion"
    STOP
  END IF

  ALLOCATE(layers_in_layer(x_extent, y_extent, z_extent), STAT=alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not allocate space for wet_proportion"
    STOP
  END IF

  ALLOCATE(layer_mass(x_extent, y_extent, z_extent, 3), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_mass"
    STOP
  END IF

  ALLOCATE(layer_attributes(x_extent, y_extent, z_extent, 3), &
                                                             STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_attributes"
    STOP
  END IF

  ALLOCATE(transmissivity(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for transmissivity"
    STOP
  END IF

  ALLOCATE(layer_storage(x_extent, y_extent, z_extent, 2), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*, *) "Model could not allocate space for layer_storage"
    STOP
  END IF

  ALLOCATE(no_layers(x_extent,y_extent), STAT = alloc_error)
  IF (alloc_error /= 0) THEN
    WRITE (*,*) "Model could not alocate space for no_layers"
    STOP
  END IF

!-------------------------------------------------------------------------------
! Section 4.0 Initialisation
!-------------------------------------------------------------------------------

  !Initialise all single-valued variables and arrays.

  !Single-valued variables.
  year_counter = 1
  sub_annual_counter = 0
  annual_tsteps = 1.0
  net_rain_tstep = 0 
  timestep = t_min / REAL(525600)
  t_step_sum_netr = 0.0
  t_step_sum_su = 0.0
  mean_t_step = 0.0
  sub_year_counter = 0
  sub_year_counter_su = 0
  output_counter = 0
  oxic_decay = 0.0
  anoxic_decay = 0.0
  x = 1
  y = 1
  z = 1
  new_net_rain = 1
  weight_z = 0
  weight_z_minus = 0

  !Arrays.
  no_layers = 0
  layer_mass = 0.0
  wet_proportion = 0.0
  col_wt_depth = 0.0
  col_wt_sum = 0.0
  col_wt_sum_su = 0.0
  col_wt_depth_av = 0.0
  col_wt_depth_av_su = 0.0
  col_mass_per_area = 0.0
  layer_attributes = 0.0
  water_change = 1.0
  wk_mean = 1.0
  transmissivity = 1.0
  layer_storage = 1.0
  water_table = 0.0
  base_altitude = 0.0
  layer_age = 0.0
  layers_in_layer = 0
  ground_water_flux = 0.0

  !Read data from files to rank-two arrays for active columns
  DO x = 1, x_extent
    DO y = 1, y_extent
      READ (040, *) activation_status(x, y)
      READ (050, *) base_altitude(x, y)
      READ (054, *) ground_water_flux(x,y)
    END DO
  END DO

  !Assign three layers to active columns - mineral, z1 and pond
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
      no_layers(x,y) = 3
      END IF
    END DO
  END DO

  !Initialise pond properties. K and s only set for layer 1 for
  !"on" and "diri" columns because additional layers added each year.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
      !Pond properties
      layer_attributes(x, y, no_layers(x,y), 1) = pond_depth !thickness
      layer_attributes(x, y, no_layers(x,y), 2) = 1000 !K
      layer_attributes(x, y, no_layers(x,y), 3) = 1.0 !s
      END IF
    END DO
  END DO

  ! Initialise the properties for a layer of mineral soil.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        !Specify layer thickness
        layer_attributes(x, y, 1, 1) = mineral
        !Calculate layer 1 k
        layer_attributes(x, y, 1, 2) = 3153.6 !Layer K = 10-4cm s-1
        !Layer porosity
        layer_attributes(x, y, 1, 3) =  0.3
        !Column elevation for layer 1 = layer thickness of layer 1
        transmissivity(x, y, 1, 1) = layer_attributes(x, y, 1, 1)
        !Calculate layer 1 transmissivity = layer 1 thickness * layer 1 k
        transmissivity(x, y, 1, 2) =  transmissivity(x, y, 1, 1) &
                                       * layer_attributes(x, y, 1, 2)
        !Set wet proportion (saturated)
        wet_proportion(x, y, 1) = 1.0
        !Initial column properties
        !Assign mass per area = layer 1 mass
        !col_mass_per_area(x, y) = layer_mass(x, y, 1, 1)
        !Set the water table to the top of layer 1 for active columns
        water_table(x, y) = transmissivity(x, y, 1, 1)
      END IF
    END DO
  END DO

  ! Initialise properties of first layer of peat profile for "on" columns.
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        !Calculations and assignments for columns with either active status or
        !Calculate layer 2 mass
        ! Equation based on Belyea and Clymo (2001) - see main model loop.
        ! Assumed that water table at surface and temp. at base value
        layer_mass(x ,y, 2, 1)  = (9.3**2.0) * 0.0001
        !layer mass = 0.008649 g cm^2
        !Assign layer 2 initial mass to be layer 2 mass
        layer_mass(x, y, 2, 2) = layer_mass(x, y, 2, 1)
        !Set layer 2 remaining mass to be 100%
        layer_mass(x, y, 2, 3) = 1.0
        !Calculate layer thickness
        layer_attributes(x, y, 2, 1) = layer_mass(x, y, 2, 1) / density
        !Calculate layer 2 K
        !Equation below from Morris et al. (2011) - see main model loop.
        layer_attributes(x, y, 2, 2) =  k_param_a * (exp (k_param_b))
        !Layer porosity
        layer_attributes(x, y, 2, 3) =  porosity
        !Column elevation for layer 2 = layer thickness of layer 1 + layer 2
        transmissivity(x, y, 2, 1) = transmissivity(x, y, 1, 1) + &
                                       layer_attributes(x, y, 2, 1)
        !Calculate layer 2 transmissivity = layer 1 thickness * layer 1 k
        transmissivity(x, y, 2, 2) =  transmissivity(x, y, 1, 1) &
                                      + (layer_attributes(x, y , 2 ,1) &
                                       * layer_attributes(x, y, 2, 2))
        !Set wet proportion (saturated)
        wet_proportion(x, y, 2) = 1.0
        !Initial column properties
        !Assign mass per area = layer 2 mass
        col_mass_per_area(x, y) = layer_mass(x, y, 2, 1)
        !Set the water table to the top of layer 2 for active columns
        water_table(x, y) = transmissivity(x, y, 2, 1)
        !Set the age of the first peat layer (i.e. z = 2) as the oldest
        layer_age(x, y, 2) = t_extent
        layers_in_layer(x, y, 2) = 1
      END IF
    END DO
  END DO

  !Set water-table height for Dirichlet condition. x1 and x_extent are
  !"neu" or "off" cells as is y_extent
  Diri_wt: DO x = 2, (x_extent -1)
    DO y = 1, (y_extent -1)
      IF (activation_status(x, y) == "diri") THEN
         water_table(x, y) =  transmissivity(x, y + 1, 2, 1)
         !Or set a value independent of the next active column e.g.
         !water_table(x, y) = 1.0 
      END IF
    END DO
  END DO Diri_wt

  !Determine the read frequency for net rainfall (annual (1), monthly (2), 
  !weekly (3))
  !Annual
  IF (net_rain_res == 1) THEN
          net_rain_tstep = annual_tsteps
  !Monthly
  ELSE IF (net_rain_res == 2) THEN
          net_rain_tstep = annual_tsteps / REAL(12)
  !Weekly
  ELSE IF (net_rain_res == 3) THEN
          net_rain_tstep = annual_tsteps / REAL(52)
  END IF

  !Timestep management. If a variable timestep is not used, one that increases
  !linearly can be used to speed up the model. The rate of increase and min and 
  !max values are currently ascertained through trial and error.
  !Set the rate of timestep increase between t_min and t_max
  t_rate = (t_max - t_min) / (t_max_years - t_min_years)
  
!-------------------------------------------------------------------------------
! Section 5.0 Main Calculations; Time Management
!-------------------------------------------------------------------------------
  CALL cpu_time(time_start)

  !1.Commence main annual loop -------------------------------------------------
  Annual_loop: DO

    !Weather
    !If annual weather, read in net rainfall
    IF (net_rain_res == 1) THEN
        READ (020, *) rainfall
    END IF
    !Display model clock on screen in units of years.
    WRITE (*,*) "Model year ", year_counter
    
    !Read annual temperature and growing season info
    READ (051, *) gsTemperature !read from file the mean growing season temperature
    READ (052, *) growingSeasonBegin !read from file the week the growing season begins
    READ (053, *) growingSeasonEnd !read from file the week the growing season ends
    
    
    !Calculate decay rate variables for current year.
    !oxic_decay = oxic_decay_base &
    !  * Q10_oxic**((temperature - base_temp) / 10.0)
    !anoxic_decay = anoxic_decay_base &
    !  * Q10_anoxic**((temperature - base_temp) / 10.0)

      !2.Start sub annual time steps -------------------------------------------
      Sub_annual_loop: DO
      !2.1 Net rainfall/temperature and counters
      !Check if a new sub-annual time step has been signalled with new_year.
      !If so read in the next values for net rainfall
      IF (net_rain_res /= 1) THEN 
        IF (new_net_rain == 1) THEN
          !Read next net rainfall
          READ (020, *) rainfall

          !Reset time counters for next trip through sub annual loop
          t_step_sum_netr = 0.0
          t_step_sum_su = 0.0
          new_net_rain = 0
        END IF
      END IF
      
      !2.X Calculate decay rate variables for current week        
      oxic_decay = oxic_decay_base &
          * Q10_oxic**((gsTemperature - base_temp) / 10.0) !using mean growing season temperature
      anoxic_decay = anoxic_decay_base &
          * Q10_anoxic**((gsTemperature - base_temp) / 10.0) !using mean growing season temperature
      
      !Count the number of sub-annual loops (≡ annual_tsteps)
      !the count is used in the average water-table calc in the annual loop
      sub_year_counter = sub_year_counter + 1
      !Summer water tables
      IF (net_rain_res /= 1) THEN 
        IF (net_rain_res == 3) THEN
          IF (sub_annual_counter >= growingSeasonBegin .AND. sub_annual_counter <= growingSeasonEnd) THEN !"summer" is now the growing season
              sub_year_counter_su = sub_year_counter_su + 1
          END IF
        ELSE IF (net_rain_res == 2) THEN
          IF (sub_annual_counter == 6 .OR. sub_annual_counter == 7 .OR. &
             sub_annual_counter == 8) THEN
              sub_year_counter_su = sub_year_counter_su + 1
          END IF
        END IF
      END IF

      !Re-initialise column properties
      col_mass_per_area = 0.0

      !Set the timestep for this iteration to be equal to the relevant timescale  
      !and reset transmax to 0.0
      !timestep = t_min / REAL(525600)
      transmax = 0.0

      !2.2 Calculate water storage and depth averaged K
      !Calculate amount of water (expressed as a depth) which may be stored in
      !each layer within each column.
      CALL layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                        layer_storage, activation_status)

      !Calculate the depth-averaged K below the water table for each column
      CALL wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)

      !2.3 Set maximum timesteps depending on elapsed model years
      !Reduce the calculated timestep for use with both hydrological and
      !decomposition routines
      !timestep = (0.5 * ((0.5 * spatial_step) ** 2) * porosity / transmax)

      !Reduce the time step if necessary
      !timestep = timestep * t_step_multi

      !Reduce the calculated timestep for use with both hydrological and
      !decomposition routines. Timestep set to t_min for first t_min_years years
      !and then to a maxumum of t_max minutes thereafter
      ! IF (year_counter > t_min_years) THEN
      !   !If the calculated timestep is greater than t_max mins, then set it to 
      !   !t_max
      !   IF (timestep > t_max / annual_tsteps) THEN
      !       timestep = t_max / annual_tsteps
      !   END IF
      ! !If the simulation has run for less than t_min_years, set the time step to 
      ! !t_min
      ! ELSE IF(year_counter <= t_min_years) THEN
      !   IF (timestep >  t_min / annual_tsteps) THEN
      !       timestep =  t_min / annual_tsteps
      !   END IF
      ! END IF
  
      !Calculate the time step
      IF (year_counter > t_min_years + t_max_years) THEN
          timestep = t_max / REAL(525600)
        ELSE IF (year_counter <=  t_min_years) THEN
          timestep = t_min / REAL(525600)
        ELSE
          timestep = (((year_counter - t_min_years) * t_rate) + t_min) / & 
                        REAL(525600)
      END IF

      !Keep the sum of sub-annual timesteps equal to the net rain interval time
      !step to prevent the over-accumulation of time in the model time extent
      IF (timestep + t_step_sum_netr > net_rain_tstep) THEN
        timestep = net_rain_tstep - t_step_sum_netr
        !For sub-annual weather, set the new net rain read flag
        IF (net_rain_res /= 1) THEN
          new_net_rain = 1
        END IF
      !Increment the sub annual counter
      sub_annual_counter = sub_annual_counter + 1
      END IF

      !Sum up the value of each timestep to keep track of sub-annual looping
      t_step_sum_netr = t_step_sum_netr + timestep
      !Summer water tables. 
      IF (net_rain_res /= 1) THEN 
        IF (net_rain_res == 3) THEN
           IF (sub_annual_counter >= growingSeasonBegin .AND. sub_annual_counter <= growingSeasonEnd) THEN !"summer" is now the growing season
            t_step_sum_su = t_step_sum_su + timestep
           END IF
        ELSE IF (net_rain_res == 2) THEN
          IF (sub_annual_counter == 6 .OR. sub_annual_counter == 7 .OR. &
                  sub_annual_counter == 8) THEN
              t_step_sum_su = t_step_sum_su + timestep
          END IF
        END IF
      END IF

      !2.4 Move water and update the water-table position
      !Calculate the amount of water (expressed as a depth) that moves between
      !each column. The flow law is based on the Boussinesq equation.
      CALL move_water(x_extent, y_extent, timestep, spatial_step, rainfall, &
                      base_altitude, water_change, water_table, wk_mean, &
                      activation_status, ground_water_flux)

      !Update position of the water table in each column
      CALL water_table_update(x_extent, y_extent, no_layers, water_change, &
                              water_table, layer_attributes, layer_storage, &
                              activation_status)

      !Calculate new water-table depth for each active column and use to
      !calculate mean water-table depth for year. The mean value is used to
      !calculate the new annual layer litter production.
      Water_table_depth: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            col_wt_depth(x, y) = transmissivity(x, y, (no_layers(x,y) - 1), 1) &
                               - water_table(x, y)
            col_wt_sum(x, y)  = col_wt_sum(x, y) + col_wt_depth(x, y)
            !Summer water tables 
            IF (net_rain_res /= 1) THEN
              IF (net_rain_res == 3) THEN
                IF (sub_annual_counter >= growingSeasonBegin & 
                        .AND. sub_annual_counter <= growingSeasonEnd) THEN !"summer" water table is now the the growing season water table
                  col_wt_sum_su(x, y) = col_wt_sum_su(x, y) + col_wt_depth(x, y)
                END IF
              ELSE IF (net_rain_res == 2) THEN 
                IF (sub_annual_counter == 6 .OR. sub_annual_counter == 7 .OR. &
                        sub_annual_counter == 8) THEN
                  col_wt_sum_su(x, y) = col_wt_sum_su(x, y) + col_wt_depth(x, y)
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO Water_table_depth

      !2.5 Decompose peat in each layer using info. from previous timestep.
      !Use proportionate decay model. Amount of peat lost depends on whether
      !peat above or below the water table. Scan upwards through all model
      !layers. Use new layer properties to update column properties.
      IF (sub_annual_counter >= growingSeasonBegin .AND. sub_annual_counter <= growingSeasonEnd) THEN !only perform decomposition during growing season
          Decompose: DO x = 1, x_extent
              DO y = 1, y_extent
                  IF (activation_status(x, y) == "on") THEN
                      !Decompose peat layers only
                      DO z = 2, (no_layers(x,y) - 1)
                          !Update layer mass including a recalcitrance function (see PM)
                          !applied only to oxic and anoxic decomposition
                          layer_mass(x, y, z, 1) = (layer_mass(x, y, z, 1) * &
                              wet_proportion(x, y, z) * &
                              (exp (- anoxic_decay * timestep * &
                              (layer_mass(x, y, z, 3) &
                              **anoxic_recal)))) + &
                              (layer_mass(x, y, z, 1) * &
                              (1 - wet_proportion(x, y, z)) * &
                              (exp (- oxic_decay * timestep * &
                              (layer_mass(x, y, z, 3)**oxic_recal))))
                          !Update layer remaining mass
                          layer_mass(x, y, z, 3) = layer_mass(x, y, z, 1) / &
                              layer_mass(x, y, z, 2)
                          !Layer thickness
                          layer_attributes(x, y, z, 1) = layer_mass(x, y, z, 1) / density
                          !Current layer mass for addition to mass per area array
                          col_mass_per_area(x, y) = col_mass_per_area(x, y) +&
                              layer_mass(x, y, z, 1)
                      END DO
                  END IF
              END DO
          END DO Decompose ! End of decomposition loop

          !Recalculate layer k
          !For peat layers only
          Layer_k: DO x = 1, x_extent
              DO y = 1, y_extent
                  IF (activation_status(x, y) == "on") THEN
                      DO z = 2, (no_layers(x,y) - 1)
                          layer_attributes(x, y, z, 2) = k_param_a * (exp (k_param_b &
                              * layer_mass(x, y, z, 3)))
                      END DO
                  END IF
              END DO
          END DO Layer_k

          !Calculate transmissivity profile for each column
          CALL trans_height(x_extent, y_extent, no_layers, layer_attributes, &
              transmissivity, activation_status)

          !Recalculate layer wet proportion
          Layer_wet_prop: DO x = 1, x_extent
              DO y = 1, y_extent
                  IF (activation_status(x, y) == "on") THEN
                      DO z = 1, (no_layers(x,y) - 1)
                          !Wet proportion
                          IF (z == 1) THEN
                              IF (water_table(x, y) >= transmissivity(x, y, z, 1)) THEN
                                  wet_proportion(x, y, 1) = 1.0
                              ELSE
                                  wet_proportion(x, y, 1) =  water_table(x, y) / &
                                      layer_attributes(x, y, 1, 1)
                              END IF
                          ELSE IF (water_table(x, y) >= transmissivity(x, y, z, 1)) THEN
                              wet_proportion(x, y, z) = 1.0
                          ELSE IF (water_table(x, y) <= transmissivity(x, y, z - 1, 1)) &
                              THEN
                              wet_proportion(x, y, z) = 0.0
                          ELSE
                              wet_proportion(x, y, z) = (water_table(x, y) -&
                                  transmissivity(x, y, z - 1, 1)) / layer_attributes(x, y, z, 1)
                          END IF
                      END DO
                  END IF
              END DO
          END DO Layer_wet_prop
      END IF
      
      !2.6 Check if sub annual total has been reached
      IF (net_rain_res == 3 .AND. sub_annual_counter == 52) THEN
        EXIT
      ELSE IF (net_rain_res == 2 .AND. sub_annual_counter == 12) THEN 
        EXIT
      ELSE IF (net_rain_res == 1 .AND. sub_annual_counter == 1) THEN 
        !Reset time step sum for annual weather 
        t_step_sum_netr = 0
        EXIT
      END IF
 
      !End sub-annual loop
      END DO Sub_annual_loop 
      !-------------------------------------------------------------------------


    !3.Layer Lumping
    !3.1 Lumping routine
    IF (year_counter > layer_exclusion + 3) THEN
      DO x = 1, x_extent !x extent do
        DO y = 1, y_extent ! y extent do
          lump_index = 2
          IF (activation_status(x, y) == "on") THEN
            DO z = 3, (no_layers(x,y) - 1) !z extent do
              IF (z > (no_layers(x,y) - layer_exclusion)) THEN
                lump_index = lump_index + 1
                layer_attributes(x,y,lump_index,1) = layer_attributes(x,y,z,1)
                layer_mass(x,y,lump_index,1) = layer_mass(x,y,z,1)
                layer_mass(x,y,lump_index,2) = layer_mass(x,y,z,2)
                layer_mass(x,y,lump_index,3) = layer_mass(x,y,z,3)
                layer_age(x,y,lump_index) = layer_age(x,y,z)
                layers_in_layer(x,y,lump_index) = layers_in_layer(x,y,z)
              ELSE IF ( ( layer_attributes(x,y,z,1) < threshold ) .AND. &
                (layer_attributes(x,y,lump_index,1) < threshold ) ) THEN
                weight_z = layer_mass(x,y,z,1) / (layer_mass(x,y,z,1) + &
                        layer_mass(x,y,lump_index,1) )
                weight_z_minus = layer_mass(x,y,lump_index,1) / &
                        (layer_mass(x,y,z,1) + layer_mass(x,y,lump_index,1))
                layer_attributes(x,y,lump_index,1) = &
                        layer_attributes(x,y,z,1) + &
                        layer_attributes(x,y,lump_index,1)
                layer_mass(x,y,lump_index,1) = layer_mass(x,y,z,1) + &
                        layer_mass(x,y,lump_index,1)
                layer_mass(x,y,lump_index,2) = layer_mass(x,y,z,2) + &
                        layer_mass(x,y,lump_index,2)
                layer_mass(x,y,lump_index,3) = layer_mass(x,y,lump_index,1) / &
                          layer_mass(x,y,lump_index,2)
                layer_age(x,y,lump_index) = ( ( weight_z * layer_age(x,y,z)) + &
                (weight_z_minus * layer_age(x,y,lump_index)))
                layers_in_layer(x,y,lump_index) = layers_in_layer(x,y,z) + &
                        layers_in_layer(x,y,lump_index)
                !WRITE (*,*) "Layer was lumped"
              ELSE
                lump_index = lump_index + 1
                layer_attributes(x,y,lump_index,1) = layer_attributes(x,y,z,1)
                layer_mass(x,y,lump_index,1) = layer_mass(x,y,z,1)
                layer_mass(x,y,lump_index,2) = layer_mass(x,y,z,2)
                layer_mass(x,y,lump_index,3) = layer_mass(x,y,z,3)
                layer_age(x,y,lump_index) = layer_age(x,y,z)
                layers_in_layer(x,y,lump_index) = layers_in_layer(x,y,z)
              END IF
            END DO ! end z extent do
            lump_index = lump_index + 1
            !Move ponding layer here
            layer_attributes(x,y,lump_index,1) = & 
                    layer_attributes(x,y,no_layers(x,y),1)
            layer_attributes(x,y,lump_index,2) = & 
                    layer_attributes(x,y,no_layers(x,y),2)
            layer_attributes(x,y,lump_index,3) = & 
                    layer_attributes(x,y,no_layers(x,y),3)
            !Blank Ponding layer
            wet_proportion(x,y,lump_index) = 0.0
            layer_age(x,y,lump_index) = 0.0
            layers_in_layer(x,y,lump_index) = 0.0
            layer_mass(x,y,lump_index,1) = 0.0
            layer_mass(x,y,lump_index,2) = 0.0
            layer_mass(x,y,lump_index,3) = 0.0
            layer_storage(x,y,lump_index,1) = 0.0
            layer_storage(x,y,lump_index,2) = 0.0
            transmissivity(x,y,lump_index,1) = 0.0
            transmissivity(x,y,lump_index,2) = 0.0
            !Array blanking
            wet_proportion(x,y,z) = 0.0
            layer_age(x,y,z) = 0.0
            layers_in_layer(x,y,z) = 0.0
            layer_mass(x,y,z,1) = 0.0
            layer_mass(x,y,z,2) = 0.0
            layer_mass(x,y,z,3) = 0.0
            layer_attributes(x,y,z,1) = 0.0
            layer_attributes(x,y,z,2) = 0.0
            layer_attributes(x,y,z,3) = 0.0
            layer_storage(x,y,z,1) = 0.0
            layer_storage(x,y,z,2) = 0.0
            transmissivity(x,y,z,1) = 0.0
            transmissivity(x,y,z,2) = 0.0
            no_layers(x,y) = lump_index
          END IF !end activation status check
        END DO! end y extent do
      END DO !end x extent do

      !3.2 Re-calculate layer K, water storage, depth averaged K and 
      !transmissivity. Calculate amount of water (expressed as a depth) which 
      !may be stored in each layer within each column.

      !Recalculate layer K
      !For peat layers only
      Layer_k_lump: DO x = 1, x_extent
        DO y = 1, y_extent
          IF (activation_status(x, y) == "on") THEN
            DO z = 2, (no_layers(x,y) - 1)
            layer_attributes(x, y, z, 2) = k_param_a * (exp (k_param_b &
                                           * layer_mass(x, y, z, 3)))
            END DO
          END IF
        END DO
      END DO Layer_k_lump
  
      CALL layer_water_depth(x_extent, y_extent, no_layers, layer_attributes, &
                            layer_storage, activation_status)
      CALL wat_k_mean(x_extent, y_extent, no_layers, transmax, water_table, &
                      wk_mean, layer_attributes, transmissivity, &
                      activation_status)
      CALL trans_height(x_extent, y_extent, no_layers, layer_attributes, &
                      transmissivity, activation_status)
    END IF

    !4.Write the number of layers in a layer
    DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
        WRITE (180,*) no_layers(x,y)
        END IF
      END DO
    END DO

    !5.Calculate the average annual and summer water-table depths to be used for
    !the next layer litter calculation
    Water_table_depth_av: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
           col_wt_depth_av(x, y) = col_wt_sum(x, y) / (REAL(sub_year_counter))
           IF (net_rain_res /= 1) THEN
             col_wt_depth_av_su(x, y) = col_wt_sum_su(x, y) / &
                                      (REAL(sub_year_counter_su))
           END IF
        END IF
      END DO
    END DO Water_table_depth_av

    !6.Calculate the mean timestep in minutes and reset counters
    mean_t_step = (1.0 / REAL(sub_year_counter)) * 525600
    WRITE(130, '(20F20.8)') mean_t_step
    !Reset mean timestep for the next main loop
    mean_t_step = 0.0
    !Reset sub annual counter for new year 
    sub_annual_counter = 0

    !7.Check for model write out
    output_counter = output_counter + 1
    IF (output_counter == output_interval) THEN
      !Write results to file
      DO x = 1, x_extent
        DO y = 1, y_extent
          IF(activation_status (x,y) == "off" &
            .OR. activation_status (x,y) == "diri" &
            .OR. activation_status (x,y) == "neu") THEN
            !Column height
            WRITE(070, '(20F20.8)') -999.0
            !Water-table height
            WRITE(080, '(20F20.8)') -999.0
            !Transmissivity
            WRITE(090, '(20F20.8)') -999.0
            !Column mass per area
            WRITE(100, '(20F20.8)') -999.0
            !Column mean water-table depth
            WRITE(110, '(20F20.8)') -999.0
            !Summer mean column water-table depth
            WRITE(160, '(20F20.8)') -999.0
          ELSE
            WRITE(070, '(20F20.8)') &
            transmissivity(x, y, (no_layers(x,y) - 1), 1)
            WRITE(080, '(20F20.8)') water_table(x, y)
            WRITE(090, '(20F20.8)') transmissivity(x, y, no_layers(x,y) - 1, 2)
            WRITE(100, '(20F20.8)') col_mass_per_area(x, y)
            WRITE(110, '(20F20.8)') col_wt_depth_av(x, y)
            WRITE(160, '(20F20.8)') col_wt_depth_av_su(x, y)
          END IF
        END DO
      END DO
    END IF

    !8.Annual counting
    !Update counter for annual update
    year_counter = year_counter + 1 ! Counts annual timesteps

    !Exit main annual loop if year counter exceeds model timespan
    IF (year_counter > t_extent) EXIT

    !Reset column water-table sum for next annual loop
    col_wt_sum = 0.0
    col_wt_sum_su = 0.0
    !Reset sub-year counting for the next annual loop
    sub_year_counter = 0
    sub_year_counter_su = 0

    !9.Add new layer for each additional year and initialise new properties
    !current layer counter +1 for pond
    DO x = 1, x_extent
      DO y = 1, y_extent
        no_layers(x,y) = no_layers(x,y) + 1
      END DO
    END DO

    !10.Initialise above-ground storage layer properties for new ponding layer
    !for active columns only
    Pond_props: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
         !Initialise the properties of the pond
          layer_attributes(x, y, no_layers(x,y), 1) = pond_depth
          layer_attributes(x, y, no_layers(x,y), 2) = 0.0 !k
          layer_attributes(x, y, no_layers(x,y), 3) = 1.0 !s
        END IF
      END DO
    END DO Pond_props

    !11.Calculate mass of new annual litter layer using modification of
    !empirical quartic function in Belyea and Clymo (2001) Proceedings
    !Royal Society London B, 268, 1315-1321. The modification takes
    !account of temperature effects.

    !Is average water table > 66.8 cm below surface? If so, assign v. low
    !productivity to new annual layer (= mass addition because timesteps
    !for new layer addition are annual).
    New_layer_props: DO x = 1, x_extent
      DO y = 1, y_extent
        IF (activation_status(x, y) == "on") THEN
          !Layer age of most recent layer
          layer_age(x, y, no_layers(x, y) -1) = t_extent - year_counter + 1 
          !layers in layer
          layers_in_layer(x, y, no_layers(x, y) -1) = 1
          !Layer mass of new top layer
          IF (col_wt_depth_av_su(x, y) > 66.8) THEN !using average summer (growing season) water tables
            layer_mass(x, y, (no_layers(x,y) - 1), 1) = 0.0000001
          ELSE
            !this is the ponding depth (cm) when summer water tables are -ve.
            IF (col_wt_depth_av_su(x, y) < 0) THEN
              col_wt_depth_av_su(x, y) = 0
            END IF
            layer_mass(x, y, (no_layers(x,y) - 1), 1) =& !using average summer (growing season) water tables
                     ((9.3 + (1.33 * col_wt_depth_av_su(x, y)) - (0.022&
                     * (col_wt_depth_av_su(x, y)**2.0)))**2.0) * 0.0001&
                     *((0.1575 * gsTemperature) + 0.009132) * ((growingSeasonEnd - growingSeasonBegin) / 52.0)!changed temperature to mean growing season temperature and production reduced according to length of growing season
          END IF
          WRITE(230, '(20F20.8)') layer_mass(x, y, (no_layers(x,y) - 1), 1)  !save annual production to file
          !Layer initial mass of new top layer = layer mass of top layer
          layer_mass(x, y, (no_layers(x,y) - 1), 2) = &
                layer_mass(x, y, (no_layers(x,y) - 1), 1)
          !Layer remaining mass of new top layer = 100%
          layer_mass(x, y, (no_layers(x,y) - 1), 3) = 1.0
          !Layer thickness of new top layer
          layer_attributes(x, y, (no_layers(x,y) - 1), 1) = &
            layer_mass(x, y, (no_layers(x,y) - 1), 1) / density
          !Layer drainable porosity
          layer_attributes(x, y, (no_layers(x,y) - 1), 3) = porosity
          !Layer elevation of new top layer =
          transmissivity(x, y, (no_layers(x,y) - 1), 1) = &
          !layer elevation below new top layer + thickness of new layer
          transmissivity(x, y, (no_layers(x,y) - 1) - 1, 1) + &
              layer_attributes(x, y, (no_layers(x,y) - 1), 1)
          !Wet proportion of new top layer
          IF (water_table(x, y) >= &
                  transmissivity(x, y, (no_layers(x, y) - 1), 1)) THEN
              wet_proportion(x, y, (no_layers(x,y) - 1)) = 1.0
          ELSE IF (water_table(x, y) <= &
              transmissivity(x, y, (no_layers(x,y) - 1), 1)) THEN
                wet_proportion(x, y, z) = 0.0
          ELSE
            wet_proportion(x, y, (no_layers(x,y) -1)) = (water_table(x, y) - &
              transmissivity(x, y, (no_layers(x,y) - 1), 1)) / & 
              layer_attributes(x, y, (no_layers(x,y) - 1), 1)
          END IF
          !Set layer Ksat if new layer contains water
          IF (wet_proportion(x, y, (no_layers(x,y) - 1)) > 0.0) THEN
              layer_attributes(x, y, (no_layers(x,y) - 1), 2) =  &
                k_param_a * (exp (k_param_b))
          !Calculate transmissivity with new layer
          transmissivity(x, y, (no_layers(x,y) - 1), 2) = &
              transmissivity(x, y, (no_layers(x,y) - 1) - 1, 2) &
              + (layer_attributes(x, y, (no_layers(x,y) -1), 1) &
              * layer_attributes(x, y, (no_layers(x,y) - 1), 2))
          END IF
          col_mass_per_area(x, y) = col_mass_per_area(x, y) +&
            layer_mass(x, y, (no_layers(x,y) - 1), 1)
        END IF
      END DO
    END DO New_layer_props

    !12. Annual loop management
    !Re-set output counter
    output_counter = 0

  !End of main annual loop
  END DO Annual_loop 
  !-----------------------------------------------------------------------------

  !13.Write array outputs
  DO x = 1, x_extent
    DO y = 1, y_extent
      IF (activation_status(x, y) == "on") THEN
        DO z = 2, (no_layers(x,y) -1)
          WRITE(060, '(20F20.8)') layer_mass(x, y, z, 3)     !mass remaining
          WRITE(200, '(20F20.8)') layer_mass(x, y, z, 2)     !initial mass
          WRITE(210, '(20F20.8)') layer_attributes(x, y, z, 1) !layer thickness
          WRITE(120, '(20F20.8)') wet_proportion(x, y, z)    !layer wet prop_n
          WRITE(190,*) y, ",", z, ",", layer_mass(x, y, z, 3), ",", &
            layer_attributes(x,y,z,1), ",", transmissivity(x,y,z,1)
          WRITE(220, 1100) x, y, z, layers_in_layer(x, y, z), &!number of layers
            layer_age(x, y, z), & !layer age
            layer_attributes(x, y, z, 1), & ! layer thickness
            layer_mass(x, y, z, 2),  layer_mass(x, y, z, 3) !layer mass
        END DO
      END IF
    END DO
  END DO

  !14.Write the model elapsed time
  CALL cpu_time(time_finish)
  time_elapsed = time_finish - time_start
  WRITE(170, *) time_elapsed

END PROGRAM DIGIBOG_3D
