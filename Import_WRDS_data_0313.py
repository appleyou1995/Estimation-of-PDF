from sqlalchemy import text
import wrds

import os
import pandas as pd
import numpy  as np


# %%  Connect to WRDS & Get the properties of database [ Option ]

conn = wrds.Connection(wrds_username='irisyu')

libraries          = conn.list_libraries()
tables_optionm_all = conn.list_tables(library='optionm_all')
col_headers_option = conn.describe_table(library='optionm_all', table='opprcd2013')

query_option = text("""
                    SELECT 
                        secid, date, exdate, cp_flag, strike_price,
                        best_bid, best_offer, open_interest, volume, impl_volatility, 
                        delta, gamma, theta, vega, optionid
                    FROM 
                        optionm_all.opprcd2013
                    WHERE            
                        secid = '108105'
                    AND
                        date = '2013-03-13'
                    """)

df_20130313_raw = conn.raw_sql(query_option)

print(df_20130313_raw.head())

conn.close()


# %%  Connect to WRDS & Get the properties of database [ Index Dividend ]

conn = wrds.Connection(wrds_username='irisyu')

tables_optionm_all

# OptionMetrics - Index Dividend Yield
col_headers_optionm = conn.describe_table(library='optionm_all', table='idxdvd')

query_optionm = text("""
                     SELECT 
                         secid, date, rate
                     FROM 
                         optionm_all.idxdvd
                     WHERE 
                         date BETWEEN '2013-03-01' AND '2013-03-31'
                     AND
                         secid = '108105'
                     """)
                   
df_optionm_div = conn.raw_sql(query_optionm)
df_optionm_div['dividend_yield'] = df_optionm_div['rate'] / 100

df_optionm_div['date'] = pd.to_datetime(df_optionm_div['date'])
df_optionm_div = df_optionm_div.set_index('date')


# import another data source
Path_Input  = 'D:/Google/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/'
IndexDivYield = pd.read_csv(os.path.join(Path_Input, 'IndexDivYield20132019.txt'), delimiter=' ', header=None)
IndexDivYield.columns = ['SECID', 'date', 'dividend_yield']
IndexDivYield['date'] = pd.to_datetime(IndexDivYield['date'].astype(str), format='%Y%m%d')
IndexDivYield_july_2013 = IndexDivYield[(IndexDivYield['date'] >= '2013-03-01') & (IndexDivYield['date'] <= '2013-03-31')]
IndexDivYield_july_2013 = IndexDivYield_july_2013.set_index('date')


df_optionm_div['check'] = IndexDivYield_july_2013['dividend_yield']

conn.close()


# %%  Connect to WRDS & Get the properties of database [ Risk-free Rate ]

conn = wrds.Connection(wrds_username='irisyu')

tables_optionm_all

# OptionMetrics - Zero Coupon Yield Curve
conn.describe_table(library='optionm_all', table='zerocd')

query_optionm_rate = text("""
                          SELECT
                              date, days, rate
                          FROM
                              optionm_all.zerocd
                          WHERE 
                              date = '2013-03-13'
                          """)
                          
df_rate = conn.raw_sql(query_optionm_rate)
df_rate['rate'] = df_rate['rate'] / 100

df_rate['date'] = pd.to_datetime(df_rate['date'])
df_rate = df_rate.set_index('date')

conn.close()
                   
                   
# %%  Setting query & Load data [ S&P 500 Index ]

conn = wrds.Connection(wrds_username='irisyu')

query_SP500 = text("""
                   SELECT caldt, spindx
                   FROM crsp.dsp500
                   WHERE caldt = '2013-03-13'
                   """)
                    
df_SP500 = conn.raw_sql(query_SP500)
spindx = df_SP500['spindx'].values[0]

conn.close()


# %%  Data preprocessing

df_20130313 = df_20130313_raw
df_20130313 = df_20130313.dropna(subset=['impl_volatility'])
df_20130313['strike_price'] = df_20130313['strike_price'] / 1000

# Add S&P500 Index level
df_20130313.insert(df_20130313.columns.get_loc('strike_price') + 1, 'spindx', spindx)

# Add Weekday
df_20130313['date'] = pd.to_datetime(df_20130313['date'])
df_20130313['exdate'] = pd.to_datetime(df_20130313['exdate'])
df_20130313['weekday'] = df_20130313['exdate'].dt.day_name()

# Correction of Expiration Date: Saturday to Friday
df_20130313.loc[df_20130313['weekday'] == 'Saturday', 'exdate'] -= pd.DateOffset(days=1)
df_20130313['weekday'] = df_20130313['exdate'].dt.day_name()

# Time-to-Maturity
df_20130313['TTM'] = (df_20130313['exdate'] - df_20130313['date']).dt.days

# Add dividend yield
df_20130313['dividend_yield'] = df_20130313['date'].map(df_optionm_div.set_index(df_optionm_div.index)['dividend_yield'])

# Add risk-free rate by one-dimensional linear interpolation
df_20130313['risk_free_rate'] = np.interp(df_20130313['TTM'], df_rate['days'], df_rate['rate'])


# %%  Filtering

# bid > 0.375
df_20130313_filter = df_20130313[df_20130313['best_bid'] > 0.375]

# Time-to-Maturity = 1 month
df_20130313_filter = df_20130313_filter[df_20130313_filter['TTM'] == 30]

# Out-of-the-Money (OTM)
df_20130313_filter['keep'] = ""

df_20130313_filter.loc[
    ((df_20130313_filter['spindx'] > df_20130313_filter['strike_price']) & (df_20130313_filter['cp_flag'] == 'P')) |
    ((df_20130313_filter['spindx'] < df_20130313_filter['strike_price']) & (df_20130313_filter['cp_flag'] == 'C')) |
    (df_20130313_filter['spindx'] == df_20130313_filter['strike_price']),
    'keep'
] = 'True'

df_20130313_filter = df_20130313_filter[df_20130313_filter['keep'] == 'True']

df_20130313_filter = df_20130313_filter.sort_values(
    by=['date', 'exdate', 'strike_price', 'spindx', 'best_bid'],
    ascending=True)

# Output
current_directory = os.getcwd()
output_file = os.path.join(current_directory, '20130313_filter.csv')
df_20130313_filter.to_csv(output_file, index=True)


# %%  Double Check

Path_Input  = '/Users/irisyu/Library/CloudStorage/GoogleDrive-jouping.yu@gmail.com/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/IndexOptions20132019_SP500/'
Path_Input  = 'D:/Google/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/IndexOptions20132019_SP500/'

OP2013 = pd.read_csv(os.path.join(Path_Input, 'OP2013_7.txt'), delimiter=' ', header=None)

OP_20130313 = OP2013[OP2013.iloc[:, 1] == 20130313]
OP_20130313.columns = ['secid', 'date', 'TTM', 'cp_flag', 'strike_price', 'forward_price', 'S0',
       'best_bid', 'best_offer', 'open_interest', 'volume', 'impl_volatility',
       'delta', 'gamma', 'theta', 'vega']

OP_20130313_sorted = OP_20130313.sort_values(by='impl_volatility')
