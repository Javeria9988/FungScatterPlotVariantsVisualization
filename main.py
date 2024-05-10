import pandas as pd
import plotly.express as px
info = pd.read_excel("E:/distance plot/final611describe.xlsx")

# Data cleaning
info = info.dropna(subset=["Gene Start", "Position", "Gene End"])

# Columns data conversion to numbers
cols = ["Gene Start", "Gene End", "Position"]
info[cols] = info[cols].apply(pd.to_numeric)

# For positive strand distance calculation
def threedistance(data):
    if data['Strand'] == '+':
        return data['Gene End'] - data['Position']
    else:
        return data['Position'] - data['Gene Start']

def fivedistance(data):
    if data['Strand'] == '+':
        return data['Position'] - data['Gene Start']
    else:
        return data['Gene End'] - data['Position']

info['Variant_3prime_Distance'] = info.apply(threedistance, axis=1)
info['Variant_5prime_Distance'] = info.apply(fivedistance, axis=1)

# For negative strand
index_negate = info['Strand'] == '-'

def negativethree(data):
    if pd.notnull(data['Gene End']) and pd.notnull(data['Position']):
        return data['Position'] - data['Gene End']
    else:
        return pd.NA

def negativefive(data):
    if pd.notnull(data['Gene Start']) and pd.notnull(data['Position']):
        return data['Gene Start'] - data['Position']
    else:
        return pd.NA

info.loc[index_negate, 'Variant_3prime_Distance'] = info.loc[index_negate].apply(negativethree, axis=1)
info.loc[index_negate, 'Variant_5prime_Distance'] = info.loc[index_negate].apply(negativefive, axis=1)

# Graph
fig = px.scatter(info, x='Variant_3prime_Distance', y='Variant_5prime_Distance', color='Chromosome', title='Af_6_1_1 Variants Distance to Next Gene',
                 labels={'Variant_3prime_Distance': "Variant 3' Distance to Next Gene", 'Variant_5prime_Distance': "Variant 5' Distance from Next Gene", 'Chromosome': 'Chromosome'})
fig.show()
