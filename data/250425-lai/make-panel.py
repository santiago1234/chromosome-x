import pandas as pd

# Load the population info and the native panel
popinfo = pd.read_csv('../250320-1tgp/data/integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
native_panel = pd.read_csv('../250320-1tgp/data/native-american.txt', sep='\t', header=None)

# Prepare native panel
native_panel = native_panel.rename(columns={0: '#Sample'})
native_panel['Panel'] = 'AMR'

# Select European and African reference populations
eur_panel = popinfo[popinfo['pop'] == 'GBR'].copy()
afr_panel = popinfo[popinfo['pop'] == 'YRI'].copy()

# Merge all panels
panel = pd.concat([eur_panel, afr_panel])
panel = panel[['sample', 'super_pop']]
panel = panel.rename(columns={'sample': '#Sample', 'super_pop': 'Panel'})
panel = pd.concat([panel, native_panel])

# Save to space-separated table
panel.to_csv('data/panel.txt', sep='\t', index=False)

# Extract query population samples
query = popinfo[popinfo['pop'].isin(['MXL', 'ASW', 'PUR', 'PEL'])]

# Save query sample list without header or index
query['sample'].to_csv('data/query_samples.txt', sep=' ', index=False, header=False)

