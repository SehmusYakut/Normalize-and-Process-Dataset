import os
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn import metrics, svm
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix, classification_report

# Read the dataset
df = pd.read_csv(r"geo_downloads\GSE12195\GSE12195.csv", index_col=0)

# Show descriptive statistics of the dataset
print(df.describe())

# Read and process the labels
df_label = pd.read_csv(r"geo_downloads\GSE12195\GSE12195_label.csv")
df_label['x'] = df_label['x'].str.split(',').str[0].str.split('_').str[0]

# Obtain unique values and assign IDs
unique_values = df_label['x'].unique()
id_mapping = {value: idx for idx, value in enumerate(unique_values)}

# Create df_change variable using IDs
df_label['x'] = df_label['x'].map(id_mapping)

# Convert labels to numerical values
label_encoder = LabelEncoder()
df_label_encoded = label_encoder.fit_transform(df_label['x'])

# Visualize label distribution
stage_counts = pd.Series(df_label_encoded).value_counts()
plt.figure(figsize=(10, 6))
bar_plot = stage_counts.plot(kind='bar')
plt.title('Cell Distributions')
plt.xlabel('Cell')
plt.ylabel('Count')
plt.xticks(rotation=45)
for index, value in enumerate(stage_counts):
    plt.text(index, value, str(value), ha='center', va='bottom')
plt.show()

#Train-Test Split
# Convert categorical data to numerical values using LabelEncoder
label_encoder = LabelEncoder()
df_encoded = df.apply(label_encoder.fit_transform)

# Convert df_label to a 1-dimensional array
df_label_1d = np.ravel(df_label)

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(df_encoded, df_label_1d, test_size=0.2, random_state=0, stratify=df_label_1d)

# Define the SVM model
clf = svm.SVC(kernel="linear")

# Train the model with training data
clf.fit(X_train, y_train)

# Make predictions for the test dataset
y_pred = clf.predict(X_test)

print("X_train:", X_train.shape, "X_test", X_test.shape)
print("y_train:", y_train.shape, "y_test", y_test.shape)

#SVM
print("Accuracy:",metrics.accuracy_score(y_test, y_pred))
print("Macro Precision:",metrics.precision_score(y_test, y_pred,average='macro'))
print("Macro Recall:",metrics.recall_score(y_test, y_pred,average='macro'))

target_names = unique_values
target_names

print(classification_report(y_test, y_pred, target_names=target_names))

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(df_encoded, df_label_encoded, test_size=0.2, random_state=0, stratify=df_label_encoded)

# Perform hyperparameter optimization for Random Forest
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 5, 10]
}
grid_search = GridSearchCV(RandomForestClassifier(random_state=42), param_grid, cv=5)
grid_search.fit(X_train, y_train)
best_params_rf = grid_search.best_params_

# Create and train the Random Forest model with the best parameters
clf_rf = RandomForestClassifier(**best_params_rf, random_state=42)
clf_rf.fit(X_train, y_train)

# Make predictions for the test dataset
y_pred_rf = clf_rf.predict(X_test)

# Calculate performance metrics for Random Forest
accuracy_rf = metrics.accuracy_score(y_test, y_pred_rf)
precision_rf = metrics.precision_score(y_test, y_pred_rf, average='macro', zero_division=0)
recall_rf = metrics.recall_score(y_test, y_pred_rf, average='macro', zero_division=0)

print("Random Forest Accuracy:", accuracy_rf)
print("Random Forest Macro Precision:", precision_rf)
print("Random Forest Macro Recall:", recall_rf)

# Print classification report for Random Forest
target_names_rf = label_encoder.inverse_transform(np.unique(df_label_encoded))
print("Random Forest Classification Report:")
print(classification_report(y_test, y_pred_rf, target_names=target_names_rf))

# Calculate confusion matrix for Random Forest
conf_matrix_rf = confusion_matrix(y_test, y_pred_rf)

# Visualize confusion matrix for Random Forest
plt.figure(figsize=(10, 6))
sns.heatmap(conf_matrix_rf, annot=True, fmt='d', cmap='Blues')
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Random Forest Confusion Matrix')
plt.show()

# Convert labels to numerical values
label_encoder = LabelEncoder()
df_label_encoded = label_encoder.fit_transform(df_label['x'])

# Split the dataset into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(df_encoded, df_label_encoded, test_size=0.2, random_state=0, stratify=df_label_encoded)

# Create and train the Logistic Regression model
model_lr = LogisticRegression(max_iter=10000)
model_lr.fit(X_train, y_train)

# Make predictions with Logistic Regression
y_pred_lr = model_lr.predict(X_test)

# Calculate performance metrics for Logistic Regression
accuracy_lr = metrics.accuracy_score(y_test, y_pred_lr)
precision_lr = metrics.precision_score(y_test, y_pred_lr, average='macro', zero_division=0)
recall_lr = metrics.recall_score(y_test, y_pred_lr, average='macro', zero_division=0)

print("Logistic Regression Accuracy:", accuracy_lr)
print("Logistic Regression Macro Precision:", precision_lr)
print("Logistic Regression Macro Recall:", recall_lr)

# Print classification report for Logistic Regression
target_names_lr = [str(i) for i in np.unique(df_label_encoded)]  # Convert numpy.int64 to strings
print("Logistic Regression Classification Report:")
print(classification_report(y_test, y_pred_lr, target_names=target_names_lr))

# Calculate confusion matrix for Logistic Regression
conf_matrix_lr = confusion_matrix(y_test, y_pred_lr)

# Visualize confusion matrix for Logistic Regression
plt.figure(figsize=(10, 6))
sns.heatmap(conf_matrix_lr, annot=True, fmt='d', cmap='Blues')
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Logistic Regression Confusion Matrix')
plt.show()
