# Data Reader

The `data_reader` module provides a unified interface for reading sequence and label data 
from different file formats.  

Currently supported file types:

- **CSV files** (`.csv`)
- **Excel files** (`.xls`, `.xlsx`)
- **TXT files** (`.txt`, tab-delimited)

---

## Data Format Requirements

Input files must contain at least two columns:

- `Seqs`: stores sequence data (a single sequence or a list of sequences)
- `Labels`: stores corresponding labels (a single label or a list of labels)

---

## Example Formats

### 1. Single sequence, single label
```text
Seqs    Labels
ATAC... 1
or
'ATAC...'    '1'
```

---

### 2. Multiple sequences, single label
```text
Seqs    Labels
['ATAC...','ATAC...']    1
```

---

### 3. Single sequence, multiple labels
```text
Seqs    Labels
ATAC...    ['1','1']
```

---

### 4. Multiple sequences, multiple labels
```text
Seqs    Labels
['ATAC...','ATAC...']    ['1','1']

```
