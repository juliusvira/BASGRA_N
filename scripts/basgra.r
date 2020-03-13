read_weather_Bioforsk <- function(y = year_start,
                                  d = doy_start,
                                  n = NDAYS,
                                  f = file_weather) {
  df_weather            <- read.table( f, header=TRUE )
  row_start             <- 1
  while( df_weather[row_start,]$YR  < y ) { row_start <- row_start+1 }
  while( df_weather[row_start,]$doy < d ) { row_start <- row_start+1 }
  df_weather_sim        <- df_weather[row_start:(row_start+n-1),]
  NMAXDAYS              <- nrow(df_weather_sim)
  NWEATHER              <- as.integer(8)
  matrix_weather        <- matrix( 0., nrow=NMAXDAYS, ncol=NWEATHER )
  matrix_weather[1:n,1] <- df_weather_sim$YR
  matrix_weather[1:n,2] <- df_weather_sim$doy
  matrix_weather[1:n,3] <- df_weather_sim$GR
  matrix_weather[1:n,4] <- df_weather_sim$T
  matrix_weather[1:n,5] <- df_weather_sim$T
  matrix_weather[1:n,6] <- exp(17.27*df_weather_sim$T/(df_weather_sim$T+239)) *
                               0.6108 * df_weather_sim$RH / 100
  matrix_weather[1:n,7] <- df_weather_sim$RAINI
  matrix_weather[1:n,8] <- df_weather_sim$WNI   
  return(matrix_weather)
}


year_start     <- as.integer(2018)
doy_start      <- as.integer(1)
years_cycle <- as.integer(1)
num_cycles = as.integer(1)
parcol       <- 1 # for Qvidja
soilcn_option <- as.integer(2)
if_weathergen <- FALSE

NDAYS          <- as.integer(600)
file_weather   <- 'weather/qvidja_2018_2019_extend.txt'
file_params    <- 'parameters/parameters_qvidja_v2.csv'
file_yasso_param <- 'yasso/param/yasso_param_climdrv_nitr'
file_yasso_init <- 'initialisation/qvidja_yasso_0kgN_i3.NModel1'
file_yasso_weather <- 'weather/qvidja_2017_2019.yasso30d.txt'

size_calendar <-  100

calendar_fert     <- matrix( 0, nrow=size_calendar, ncol=5 )
calendar_fert[1,] <- c( 2018, 197, 60*1000/ 10000, 1, 9 ) # 16 June 2018: 60 kgN/ha, organic, C:N=9
calendar_fert[2,] <- c( 2018, 236, 40*1000/ 10000, 1, 9 ) # 24 August 2018
calendar_fert[3,] <- c( 2019, 121, 100*1000/ 10000, 1, 9 ) # 1 May 2019
calendar_fert[4,] <- c( 2019, 176, 50.6*1000/ 10000, 0, 0 ) # 25 June 2019: 50.6 kgN mineral

calendar_ndep     <- matrix( 0, nrow=size_calendar, ncol=3 )
calendar_ndep[1,] <- c( 1900,   1,  2*1000/(10000*365) ) # 2 kg N ha-1 y-1 N-deposition in 1900
calendar_ndep[2,] <- c( 1980, 366,  4*1000/(10000*365) ) # 4 kg N ha-1 y-1 N-deposition in 1980
calendar_ndep[3,] <- c( 2100, 366,  4*1000/(10000*365) ) # 4 kg N ha-1 y-1 N-deposition in 2100

days_harvest      <- matrix( as.integer(-1), nrow=size_calendar, ncol=3 )
days_harvest [1,] <- as.integer(c( 2018, 163, 0 )) # 12 Jun 2018 harvest
days_harvest [2,] <- as.integer(c( 2018, 233, 1 )) # 21 Aug 2018 cutting
days_harvest [3,] <- as.integer(c( 2018, 266, 0 )) # 23 Sep 2018 harvest
days_harvest [4,] <- as.integer(c( 2019, 162, 0 )) # 11 Jun 2019 harvest
days_harvest [5,] <- as.integer(c( 2019, 232, 0 )) # 20 Aug 2019 harvest

df_params      <- read.table(file_params,header=T,sep=",",row.names=1)
params         <- df_params[,parcol]
matrix_weather    <- read_weather_Bioforsk(year_start,doy_start,NDAYS,file_weather)
yasso_param <- as.matrix(scan(file_yasso_param))
yasso_init <- as.matrix(scan(file_yasso_init))
yasso_weather <- read.table(file_yasso_weather, header=FALSE, sep=',')

yasso_weather <- yasso_weather[yasso_weather[,1] > 2017,][1:600,] # yasso_weather and matrix_weather need to both have NDAYS rows

nyears <- as.integer(1) # no cycling
size_out <- as.integer(112)
matrix_out <- matrix(0, NDAYS*nyears, size_out)
dyn.load("BASGRA.DLL")
output = .Fortran('BASGRA', params, matrix_weather, t(as.matrix(yasso_weather)), # transpose to row-major
                  yasso_init, yasso_param, calendar_fert, calendar_ndep, days_harvest,
                  soilcn_option, if_weathergen, NDAYS, nyears,
                  ncol(matrix_weather), # nweather
                  nrow(matrix_weather), # size_weather
                  as.integer(size_calendar),
                  nrow(yasso_init), 
                  nrow(yasso_param),
                  as.integer(size_out), # "NOUT"
                  matrix_out)[[19]]

write.table(output, 'basgra.test', sep=",", row.names=F)
