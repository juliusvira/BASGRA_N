import argparse
import pandas as pd
import datetime as dt
import sys

def get_clim_const(filename, sep, expand=False):
    df = pd.read_csv(filename, sep=sep)
    months = []
    # create the month column
    for doy in df.doy:
        date = dt.datetime(1991, 1, 1)+dt.timedelta(days=doy-1)
        months.append(date.month)
    df['month'] = months
    tempr_monthly = df['T'].groupby(df.month).mean()
    tempr_mean = tempr_monthly.mean()
    tempr_ampl = 0.5 * (tempr_monthly.max() - tempr_monthly.min())
    precip = df.RAINI.groupby(df.YR).sum().mean() # mean of yearly precips
    if expand:
        df = pd.DataFrame({'year':df.YR, 'doy':df.doy, 'tempr':tempr_mean, 'precip':precip, 'tempr_ampl':tempr_ampl})
        return df
    else:
        return tempr_mean, precip, tempr_ampl

def get_clim_moving(filename, sep, draw=False):
    df = pd.read_csv(filename, sep=sep)
    dates = []
    for year, doy in zip(df.YR, df.doy):
        date = dt.datetime(year, 1, 1)+dt.timedelta(days=doy-1)
        dates.append(date)

    tempr = df['T']
    tempr_rolling = tempr.rolling(30).mean()
    precip_rolling = df.RAINI.rolling(30).sum()*12
    tempr_ampl = 0
    if draw:
        df['precip_rolling'] = precip_rolling
        df['date'] = dates
        import pylab
        tempr.plot()
        tempr_rolling.plot()
        pylab.show()
        df = df.set_index('date')
        df.RAINI.plot()
        df.precip_rolling.plot()
        pylab.show()
    df = pd.DataFrame({'year':df.YR, 'doy':df.doy, 'tempr':tempr_rolling, 'precip':precip_rolling, 'tempr_ampl':0})
    return df


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('file_in')
    argparser.add_argument('--output', '-o')
    argparser.add_argument('--sep', default=' ')
    argparser.add_argument('--how', choices=['const', 'const-daily', 'rolling30d'])
    argparser.add_argument('--interactive', '-i', action='store_true')
   
    args = argparser.parse_args()
    if args.how == 'const':
        tempr_mean, precip, tempr_ampl = get_clim(args.file_in, args.sep)
        clim = '%f %f %f\n' % (tempr_mean, precip, tempr_ampl)
        def write(out):
            out.write(clim)
    elif args.how == 'const-daily':
        df = get_clim_const(args.file_in, args.sep, expand=True)
        columns = 'year doy tempr precip tempr_ampl'.split()
        def write(out):
            df[columns].to_csv(out, header=False, index=False)
    elif args.how == 'rolling30d':
        df = get_clim_moving(args.file_in, args.sep, args.interactive)
        columns = 'year doy tempr precip tempr_ampl'.split()
        def write(out):
            df[columns].to_csv(out, header=False, index=False)
    else:
        raise ValueError()
    if args.output:
        with open(args.output, 'w') as output:
            write(output)
    else:
        write(sys.stdout)



