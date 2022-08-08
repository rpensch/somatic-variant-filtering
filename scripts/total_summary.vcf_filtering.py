#!/usr/bin/env python3

import pandas as pd
import sys

summary = pd.DataFrame()
for arg in sys.argv[2:]:
    if summary.empty:
        summary = pd.read_csv(arg, sep='\t')
    else:
        summary = summary.append(pd.read_csv(arg, sep='\t'))

summary.to_csv(sys.argv[1], sep='\t', index=False)