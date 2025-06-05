import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

seed = 42

X_train_bulk_path = "./data/deconv/train/train_bulk.csv"
y_train_bulk_path = "./data/deconv/train/train_bulk_trueprops.csv"

X_test_bulk_path = "./data/deconv/test/bulkified_data.csv"

desired_order = [
    "T",
    "Endothelial",
    "Fibroblast",
    "Plasmablast",
    "B",
    "Myofibroblast",
    "NK",
    "Myeloid",
    "Mast",
]


def load_data():
    # Load and prepare train data
    X_train_bulk = pd.read_csv(X_train_bulk_path).set_index("Unnamed: 0").T
    y_train_bulk = pd.read_csv(y_train_bulk_path).set_index("highLevelType").T

    # Align train data by samples
    X_train_bulk, y_train_bulk = X_train_bulk.align(y_train_bulk, axis=0)

    # Load and prepare test data
    X_test_bulk = pd.read_csv(X_test_bulk_path).set_index("Unnamed: 0").T

    return X_train_bulk, y_train_bulk, X_test_bulk


def train_model(X_train: pd.DataFrame, y_train: pd.DataFrame):
    X_train_cp = X_train
    y_train_cp = y_train

    # Split train data into training and validation sets
    X_train, X_val, y_train, y_val = train_test_split(
        X_train, y_train, test_size=0.2, random_state=seed
    )

    # Train regression model
    model = RandomForestRegressor(random_state=seed)
    model.fit(X_train, y_train)

    # Predict on the validation set
    y_valid_pred = model.predict(X_val)

    # Compute RMSE for each cell type & average RMSE (on validation data)
    rmse_per_cell_type = np.sqrt(
        mean_squared_error(y_val, y_valid_pred, multioutput="raw_values")
    )
    average_rmse = np.mean(rmse_per_cell_type)

    print("RMSE per cell type (validation data):", rmse_per_cell_type)
    print("Average RMSE (validation data):", average_rmse)

    model.fit(X_train_cp, y_train_cp)

    return model


def pred_model(model, y_train, X_test: pd.DataFrame):
    # final prediction on test data
    y_test_pred = model.predict(X_test)

    # Get correct cell type and bulk sample names
    initial_order = y_train.columns.tolist()  # Original cell type names
    bulk_samples = X_test.index.tolist()  # List of bulk sample names

    # Create the DataFrame with the original order
    pred_props = pd.DataFrame(
        y_test_pred.T,  # Transpose to have cell types as rows and bulk samples as columns
        columns=bulk_samples,  # Bulk sample names as columns
        index=initial_order,  # Original cell type order as rows
    )

    # Reorder rows to match the desired order of cell types
    pred_props = pred_props.reindex(desired_order)

    # Add an unnamed index column starting from 0
    pred_props.reset_index(inplace=True)
    pred_props.index.name = ""  # Ensure the index column has no name

    return pred_props


def export(pred_props: pd.DataFrame, path) -> None:
    pred_props.to_csv(path, index=True, index_label="")
    print(f"Submission file saved to: {path}")


def pipeline(dir):
    path = os.path.join(dir, "pred_props.csv")
    if os.path.exists(path):
        return
    X_train, y_train, X_test = load_data()
    model = train_model(X_train, y_train)
    pred_props = pred_model(model, y_train, X_test)
    export(pred_props, path)
