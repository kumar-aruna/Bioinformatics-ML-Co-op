### Bioinformatics-ML-Co-op

## Test for Bioinformatics and ML Co-Op:
## Task 01: Gene Annotation Parsing: Downloading and parsing gene annotation from the GENCODE database using R.
### Step 01: Step 01: Install the required R packages:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "S4Vectors"))
library(GenomicFeatures)
library(S4Vectors)
```
### Step 02: Downloading gene annotation from the GENCODE database.
```
# Timeout to ensure large files have time to download:
options(timeout = max(1000, getOption("timeout")))

# Function to download the gtf file from given url
download_gtf <- function(url, destination) {
  # Parameters:
  # url:Type:String, The web address from where the GTF file should be downloaded. 
  # destination: The local path where the downloaded GTF file should be saved. 
  
  # The function uses download.file() from R's base package, which downloads files over the internet.
  # mode: is set to "wb" to write the file in binary mode, ensuring the file's integrity.
  download.file(url = url, destfile = destination, mode = "wb")
}
  
# Example usage with the corrected URL (removed the leading space)
gencode_url <- "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.primary_assembly.annotation.gtf.gz"
destination_path <- "/Users/aruna/gencode.v45.annotation.gtf.gz"

# Call the corrected function
download_gtf(gencode_url, destination_path)
```

### Task 03: Create a TxDb Object
```
txdb <- makeTxDbFromGFF(file = gtf_path, format = "gtf")
```
### Task 04: Build an S4 Object Containing All Genes and Transcripts
```
library(GenomicRanges)
genes <- genes(txdb)
transcripts <- transcriptsBy(txdb, by = "gene")

gene_transcript_association <- new("list", genes = genes, transcripts = transcripts)
print(gene_transcript_association)
```
### Task 05:  Compute Statistics and Make a Histogram
```
setClass("GeneTranscriptAssociation", 
         slots = c(genes = "GRanges", transcripts = "CompressedList"))

# Instantiate the S4 object
gene_transcript_association <- new("GeneTranscriptAssociation", genes = genes, transcripts = transcripts_by_gene)

# Compute transcript counts per gene using the 'lengths' function
transcript_counts <- lengths(gene_transcript_association@transcripts)

# Now calculate the mean, minimum, and maximum
mean_transcripts <- mean(transcript_counts)
min_transcripts <- min(transcript_counts)
max_transcripts <- max(transcript_counts)

# Display the calculated values
cat("Mean number of transcripts per gene: ", mean_transcripts, "\n")
cat("Minimum number of transcripts per gene: ", min_transcripts, "\n")
cat("Maximum number of transcripts per gene: ", max_transcripts, "\n")

# Generate and display a histogram of the number of transcripts per gene
hist(transcript_counts, 
     breaks=seq(from=min(transcript_counts), to=max(transcript_counts), by=1), 
     main="Histogram of Transcripts Per Gene", 
     xlab="Number of Transcripts", 
     ylab="Frequency", 
     border="black", 
     col="skyblue")

```
### Task 06: # Save the S4 object for future use

```
saveRDS(object = gene_transcript_association, file = "~/gene_transcript_association.rds")
```


## Task 02: Boosted Decision Tree Modeling: Building a boosted decision tree for a regression task on a provided dataset using xgboost in Python.
### Step 01: Install necessary libraries::
```
import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error
from math import exp, log


```
### Step 02: Download and load the dataset
```
data_url = "http://129.10.224.71/~apaul/data/tests/dataset.csv.gz"
df = pd.read_csv(data_url, compression='gzip')

```
### Step 03: Data Preprocessing
```
def x_scale(x, p=7.5):
    return 1/p * np.log(1 + x * (np.exp(p) - 1))

def y_scale(y):
    return np.log(1 + y) if y >= 0 else -np.log(1 - y)

# Apply transformations
df['x1'] = df['x1'].apply(x_scale)
df['y1'] = df['y1'].apply(y_scale)  # Assuming you choose y1 for regression

```
### Step 4: Splitting the Data
```
X = df[['x1', 'x2', 'x3', 'x4']]
y = df['y1']  # Assuming y1 is the target

X_train, X_temp, y_train, y_temp = train_test_split(X, y, test_size=0.4, random_state=42)
X_val, X_test, y_val, y_test = train_test_split(X_temp, y_temp, test_size=0.5, random_state=42)
```

### Step 5: Splitting the Data
```
model = xgb.XGBRegressor(objective='reg:squarederror', learning_rate=0.1, max_depth=5, n_estimators=100, subsample=0.8, colsample_bytree=0.8)
model.fit(X_train, y_train, early_stopping_rounds=10, eval_set=[(X_val, y_val)], verbose=True)

```
### Step 6: Model Evaluation
```
predictions = model.predict(X_test)
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
print(f"Test RMSE: {rmse}")

```









### Author: Aruna Kumar

  

