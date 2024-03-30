### Bioinformatics-ML-Co-op

## Test for Bioinformatics and ML Co-Op:
## Task 01: Gene Annotation Parsing: Downloading and parsing gene annotation from the GENCODE database using R.
### Step 01: Install the required R packages:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GenomicFeatures", "S4Vectors"))
library(GenomicFeatures)
library(S4Vectors)

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


## Task 02: Boosted Decision Tree Modeling: Building a boosted decision tree for a regression task on a provided dataset using xgboost in Python.






### Authors: Aruna Kumar

  

