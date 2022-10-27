import datetime as dt
import pandas as pd
import matplotlib as mp

# Clean date and time of flights


def set_datetime(df, tz, unique=False) -> object:
    if unique:
        date_time = pd.unique(df.Date)
    else:
        date_time = df.Date

    date_time = pd.to_datetime(date_time)
    date_time = pd.DatetimeIndex(date_time).tz_convert(tz)
    dates = date_time.strftime('%y-%m-%d')
    times = date_time.strftime('%H:%M:%S')

    return dates, times


def datetime2num(dates, times):
    x = [dt.datetime.strptime(d, '%y-%m-%d') for d in dates]
    y = [dt.datetime.strptime(t, '%H:%M:%S') for t in times]
    y = mp.dates.datestr2num(times)

    return x, y
