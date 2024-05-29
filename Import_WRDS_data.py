from sqlalchemy import text
import wrds


# %%  Connect to WRDS

conn = wrds.Connection(wrds_username='irisyu')


# %%  Get the properties of database

libraries = conn.list_libraries()

tables = conn.list_tables(library='optionm_all')

col_headers = conn.describe_table(library='optionm_all', table='opprcd2023')
col_headers


# %%

query = text("""
        SELECT 
            secid, date, exdate, cp_flag, strike_price,
            best_bid, best_offer, impl_volatility, optionid
        FROM 
            optionm_all.opprcd1996
        WHERE            
            secid = '108105'
        AND
            date = '1996-07-17'
        """)
# %%

# 執行查詢並獲取數據
df_19960717 = conn.raw_sql(query)

# 顯示數據
print(df_19960717.head())

# 關閉連接
conn.close()


import os
import pandas as pd
Path_Input  = '/Users/irisyu/Library/CloudStorage/GoogleDrive-jouping.yu@gmail.com/我的雲端硬碟/學術｜研究與論文/論文著作/CDI Method/Data/99 姿穎學姊提供/20240417/IndexOptions19962019_SP500/'
OP1996 = pd.read_csv(os.path.join(Path_Input, 'OP1996_7.txt'), delimiter=' ', header=None)

OP_19960717 = OP1996[OP1996.iloc[:, 1] == 19960717]
