from sqlalchemy import text
import wrds

import os
import pandas as pd


# %%  Connect to WRDS & Get the properties of database

conn = wrds.Connection(wrds_username='irisyu')

libraries          = conn.list_libraries()
tables_option      = conn.list_tables(library='optionm_all')
col_headers_option = conn.describe_table(library='optionm_all', table='opprcd2023')


# %%  Setting query & Load data

query_option = text("""
                    SELECT 
                        secid, date, exdate, cp_flag, strike_price,
                        best_bid, best_offer, open_interest, volume, impl_volatility, 
                        delta, gamma, theta, vega, optionid
                    FROM 
                        optionm_all.opprcd1996
                    WHERE            
                        secid = '108105'
                    AND
                        date = '1996-07-17'
                    """)

df_19960717_raw = conn.raw_sql(query_option)

print(df_19960717_raw.head())

conn.close()


# %%  Get the properties of database

conn = wrds.Connection(wrds_username='irisyu')

tables_SP500      = conn.list_tables(library='crsp')
col_headers_SP500 = conn.describe_table(library='crsp', table='dsp500')


# %%  Setting query & Load data

query_SP500 = text("""
                   SELECT caldt, spindx
                   FROM crsp.dsp500
                   WHERE caldt = '1996-07-17'
                   """)
                    
df_SP500 = conn.raw_sql(query_SP500)
spindx = df_SP500['spindx'].values[0]

conn.close()


# %%  Data preprocessing

df_19960717 = df_19960717_raw
df_19960717 = df_19960717.dropna(subset=['impl_volatility'])
df_19960717['strike_price'] = df_19960717['strike_price'] / 1000

# Add S&P500 Index level
df_19960717.insert(df_19960717.columns.get_loc('strike_price') + 1, 'spindx', spindx)

# Time-to-Maturity
df_19960717['date'] = pd.to_datetime(df_19960717['date'])
df_19960717['exdate'] = pd.to_datetime(df_19960717['exdate'])
df_19960717['weekday'] = df_19960717['exdate'].dt.day_name()

# Correction of Expiration Date: Saturday to Friday
df_19960717.loc[df_19960717['weekday'] == 'Saturday', 'exdate'] -= pd.DateOffset(days=1)
df_19960717['weekday'] = df_19960717['exdate'].dt.day_name()
df_19960717['TTM'] = (df_19960717['exdate'] - df_19960717['date']).dt.days


# %%  Filtering

# bid > 0.375
df_19960717_filter = df_19960717[df_19960717['best_bid'] > 0.375]

# Time-to-Maturity = 1 month
df_19960717_filter = df_19960717_filter[df_19960717_filter['TTM'] == 30]

# Out-of-the-Money (OTM)
df_19960717_filter['keep'] = ""

df_19960717_filter.loc[
    ((df_19960717_filter['spindx'] > df_19960717_filter['strike_price']) & (df_19960717_filter['cp_flag'] == 'P')) |
    ((df_19960717_filter['spindx'] < df_19960717_filter['strike_price']) & (df_19960717_filter['cp_flag'] == 'C')) |
    (df_19960717_filter['spindx'] == df_19960717_filter['strike_price']),
    'keep'
] = 'True'

df_19960717_filter = df_19960717_filter[df_19960717_filter['keep'] == 'True']

df_19960717_filter = df_19960717_filter.sort_values(
    by=['date', 'exdate', 'strike_price', 'spindx', 'best_bid'],
    ascending=True)


# %%  Double Check

Path_Input  = '/Users/irisyu/Library/CloudStorage/GoogleDrive-jouping.yu@gmail.com/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/IndexOptions19962019_SP500/'
Path_Input  = 'D:/Google/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/IndexOptions19962019_SP500/'

OP1996 = pd.read_csv(os.path.join(Path_Input, 'OP1996_7.txt'), delimiter=' ', header=None)

OP_19960717 = OP1996[OP1996.iloc[:, 1] == 19960717]
OP_19960717.columns = ['secid', 'date', 'TTM', 'cp_flag', 'strike_price', 'forward_price', 'S0',
       'best_bid', 'best_offer', 'open_interest', 'volume', 'impl_volatility',
       'delta', 'gamma', 'theta', 'vega']

OP_19960717_sorted = OP_19960717.sort_values(by='impl_volatility')
