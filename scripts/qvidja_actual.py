import pandas as pd
import basgra
import argparse

year_start = 2018
doy_start = 1

file_weather = 'weather/qvidja_2017_2019.txt'
file_params = 'parameters/parameters_qvidja_v3.csv'
file_params = 'parameters/parameters_qvidja_v3.csv'

param = pd.read_csv(file_params, sep=',', index_col=0).Qvidja

file_yasso_params_soil = 'yasso/param/yasso_param_soildrv_nitr'
file_yasso_params_clim = 'yasso/param/yasso_param_climdrv_nitr'

file_init_nmodel1 = 'initialisation/qvidja_yasso_0kgN_i3.NModel1'
file_init_nmodel2 = 'initialisation/qvidja_yasso_0kgN_i3.TwoAgeTwoEffModel'

def set_fert(mdl):
    mdl.add_fert(2018, 197, 60, 'org', 9) # 16 June 2018: 60 kgN/ha, organic
    mdl.add_fert(2018, 236, 40, 'org', 9) # 24 August 2018: organic
    mdl.add_fert(2019, 121, 100, 'org', 9) # 1 May 2019: organic
    mdl.add_fert(2019, 176, 50.6, 'min')  # 25 June 2019: 50.6 kgN mineral

def set_fert_mineral(mdl):
    mdl.add_fert(2018, 197, 60, 'min')
    mdl.add_fert(2018, 236, 40, 'min')
    mdl.add_fert(2019, 121, 100, 'min')
    mdl.add_fert(2019, 176, 50.6, 'min')
    
    
def set_harv(mdl):
    mdl.add_harvest(2018, 163)
    mdl.add_harvest(2018, 233, cut_only=True)
    mdl.add_harvest(2018, 266)
    mdl.add_harvest(2019, 162)
    mdl.add_harvest(2019, 232)
    
def setup_run(mdl, filename, fert=set_fert, harv=set_harv, file_init_yasso=file_init_nmodel1, params=None):
    mdl.load_weather(file_weather)
    mdl.set_param(param)
    if params:
        mdl.load_yasso_params(params)
    elif mdl.yasso_clim_mode == 'soil':
        mdl.load_yasso_params(file_yasso_params_soil)
    else:
        mdl.load_yasso_params(file_yasso_params_clim)
        
    mdl.load_yasso_init(file_init_yasso)
    fert(mdl)
    harv(mdl)
    mdl.run(year_start, doy_start)
    mdl.store_csv(filename)

#model_const = basgra.Basgra(soilcn_option='nmodel1_soil', yasso_clim_mode='soil')
#setup_run(model_const, 'qvidja_y1soil_2018_2019.csv.gz')

# model_const = basgra.Basgra(soilcn_option='nmodel1', yasso_clim_mode='const')
# setup_run(model_const, 'qvidja_y1const_2018_2019.csv.gz')


# model_const = basgra.Basgra(soilcn_option='basgra', yasso_clim_mode='const')
# setup_run(model_const, 'qvidja_b_2018_2019.csv.gz')

model_rolling = basgra.Basgra(soilcn_option='nmodel1', yasso_clim_mode='rolling')
setup_run(model_rolling, 'qvidja_y1rolling_2018_2019.csv.gz')

model_rolling_mm = basgra.Basgra(soilcn_option='nmodel1', yasso_clim_mode='rolling')
#setup_run(model_rolling_mm, 'qvidja_y1rollingM_2018_2019.csv.gz')
setup_run(model_rolling_mm, 'qvidja_y1rollingM_2018_2019.csv.gz', params='yasso/param/yasso_param_climdrv_monthly_yearly_nitr')


# model_minfert = basgra.Basgra(soilcn_option='nmodel1', yasso_clim_mode='rolling')
# setup_run(model_minfert, 'qvidja_y1rminf_2018_2019.csv.gz', fert=set_fert_mineral)
