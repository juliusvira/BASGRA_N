import basgrafort as bg
import pandas as pd
import numpy as np
import average4yasso
outdef_def = 'conf/basgra_vars.csv'

ndep_default = np.array([[1900,   1,  2*1000/(10000*365)],
                         [1980, 366,  4*1000/(10000*365)],
                         [ 2100, 366,  4*1000/(10000*365)]])

class Basgra:
    def __init__(self, use_tempr_min_max=False, outdef=outdef_def, soilcn_option=1, yasso_clim_mode='const'):
        cnmap = {'basgra':1, 'nmodel1':2, 'nmodel2':3, 'basgra_nm1':4, 'basgra_nm2':5, 'nmodel1_soil':6}
        for num in range(1,6):
            cnmap[num] = num
        self.use_tempr_min_max = use_tempr_min_max
        self.soilcn_option = cnmap[soilcn_option]
        self.yasso_clim_mode = yasso_clim_mode
        self.load_outdef(outdef)
        self.fertlst = []
        self.harvlst = []
        self.ndep = ndep_default
        self.yasso_init = np.zeros(1)
        self.yasso_params = np.zeros(1)
        
    def load_yasso_init(self, filename):
        self.yasso_init = np.fromfile(filename, sep=' ')

    def set_yasso_init(self, values):
        self.yasso_init = values

    def load_yasso_init(self, filename):
        self.yasso_init = np.fromfile(filename, sep=' ')
        
    def load_weather(self, filename):
        self.wx = pd.read_csv(filename, sep=' ')
        #matrix_weather[1:n,6] <- exp(17.27*df_weather_sim$T/(df_weather_sim$T+239)) *
        #                       0.6108 * df_weather_sim$RH / 100
        tempr = self.wx['T']
        self.wx['VPI'] = np.exp(17.27*tempr / (tempr+239)) * 0.6108 * self.wx.RH * 0.01
        self.yassowx = self._make_yasso_weather(filename)
        
    def _make_yasso_weather(self, filename):
        if self.yasso_clim_mode == 'const':
            df = average4yasso.get_clim_const(filename, sep=' ', expand=True)
        elif self.yasso_clim_mode == 'rolling':
            df = average4yasso.get_clim_moving(filename, sep=' ')
        elif self.yasso_clim_mode == 'soil':
            df = average4yasso.get_clim_const(filename, sep=' ', expand=True)
        else:
            raise ValueError('Bad yasso_clim_mode')
        columns = 'year doy tempr precip tempr_ampl'.split()
        return df[columns]
            
    def load_param(self, filename):
        self.param = pd.read_csv(filename, sep=',', index_col=0)

    def load_outdef(self, filename):
        self.outdef = pd.read_csv(filename, index_col=0)

    def set_param(self, param):
        self.param = param

    def set_yasso_params(self, yasso_params):
        self.yasso_params = yasso_params

    def load_yasso_params(self, filename):
        self.yasso_params = np.fromfile(filename, sep=' ')
        
    def add_fert(self, year, doy, rate, fert_type, fert_cn=0):
        typeflag = {'mineral':0, 'min':0, 'organic':1, 'org':1}[fert_type]
        if typeflag == 1 and fert_cn == 0:
            raise ValueError('fert_cn needed for organic fertilizer')
        self.fertlst.append([year, doy, 0.1*rate, typeflag, fert_cn])

    def add_harvest(self, year, doy, cut_only=False):
        harv_flag = 1 if cut_only else 0
        self.harvlst.append([year, doy, harv_flag])

    def get_fert_harv_ndep(self):
        num_rec = max(len(self.harvlst), len(self.fertlst), self.ndep.shape[0])
        if num_rec > 100:
            # this might work but have to check the basgra code
            raise ValueError()
        num_rec = 100
        fert = np.zeros((num_rec, 5))
        ndep = np.zeros((num_rec, 3))
        harv = np.zeros((num_rec, 3))
        if self.fertlst:
            fert[:len(self.fertlst),:] = np.array(self.fertlst)
        if self.harvlst:
            harv[:len(self.harvlst),:] = np.array(self.harvlst)
        ndep[:self.ndep.shape[0],:] = self.ndep
        return fert, harv, ndep
            
    @staticmethod
    def from_files(filename_param, filename_weather, filenama_outdef=outdef_def, **kwargs):
        model = Basgra(**kwargs)
        model.set_param(filename_param)
        model.set_weather(filename_weather)
        model.set_outdef(filename_outdef)
        return model

    def run(self, year_start=None, doy_start=None, num_days=None):
        print(self.wx)
        if year_start is None:
            year_start = self.wx.YR[0]
        if doy_start is None:
            doy_start = self.wx.doy[0]
        if self.use_tempr_min_max:
            wxvars = 'YR doy GR TMMXI TMMNI VPI RAINI WNI'.split()
        else:
            wxvars = 'YR doy GR T T VPI RAINI WNI'.split()
        take_rows = (self.wx.YR >= year_start) | ((self.wx.YR == year_start) & self.wx.doy >= doy_start)
        if num_days is None:
            num_days = take_rows.sum()
            print('NUM DAYS', num_days, year_start)
        matrix_weather = self.wx.loc[take_rows, wxvars].values

        # note yasso takes the weather in fortran order:
        matrix_yasso_weather = self.yassowx.loc[take_rows, :].values.T
        
        print(self.wx.loc[take_rows, wxvars])

        if self.soilcn_option > 1:
            if len(self.yasso_init) < 5:
                raise ValueError('Yasso soil is requested but no initialization given')
            if len(self.yasso_params) < 2:
                raise ValueError('Yasso soil is requested but no params given')
        num_years = 1
        num_out = self.outdef.shape[0]
        fert, harvest, ndep = self.get_fert_harv_ndep()
        if not (fert.shape[0] == ndep.shape[0] == harvest.shape[0]):
            raise ValueError('fert, ndep, harvest must have equal number of rows')
        self.output = bg.basgra(self.param.values, matrix_weather, matrix_yasso_weather, self.yasso_init, self.yasso_params, 
                                fert, ndep, harvest, self.soilcn_option, False, num_days, num_years, num_out)

    def store_csv(self, filename):
        frame = pd.DataFrame(self.output, columns=self.outdef.index)
        frame.to_csv(filename)


class Schedule:
    max_entries = 100
    columns = None
    
    def __init__(self):
        self.records = []

    def add(self, *values):
        year, doy = values[0], values[1]
        if len(self.records) == self.max_entries:
            raise ValueError('Too many schedule entries')
        if doy < 1 or doy > 365:
            raise ValueError('Bad doy')
        self.check(values)
        self.records.append(values)

    def check(self, values):
        if self.columns is not None and len(values) != len(self.columns):
            raise ValueError('Wrong number of columns')

    def values(self):
        shape = self.max_entries, len(self.columns)
        matrix = np.zeros(shape)
        matrix[0:len(self.records),:] = np.array(self.records)
        return matrix
    
class Fert(Schedule):
    def check(self):
        self.super.check()
