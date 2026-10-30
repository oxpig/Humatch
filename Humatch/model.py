import os
import requests
from tensorflow import keras
from Humatch.utils import HEAVY_V_GENE_CLASSES, LIGHT_V_GENE_CLASSES, PAIRED_CLASSES


# cnn weights added to compiled env package_data
HUMATCH_CODE_DIR = os.path.dirname(os.path.abspath(__file__))
# set HUMATCH_WEIGHTS_DIR to keep trained_models/ and germline_likeness_lookup_arrays/ elsewhere
HUMATCH_WEIGHTS_DIR = os.environ.get("HUMATCH_WEIGHTS_DIR") or HUMATCH_CODE_DIR
CNN_WEIGHTS_DIR = os.path.join(HUMATCH_WEIGHTS_DIR, "trained_models")
if not os.path.exists(CNN_WEIGHTS_DIR):
    os.makedirs(CNN_WEIGHTS_DIR)
HEAVY_WEIGHTS = os.path.join(CNN_WEIGHTS_DIR, "heavy.weights.h5")
LIGHT_WEIGHTS = os.path.join(CNN_WEIGHTS_DIR, "light.weights.h5")
PAIRED_WEIGHTS = os.path.join(CNN_WEIGHTS_DIR, "paired.weights.h5")
model_urls = {
    "heavy":  "https://zenodo.org/records/13764771/files/heavy.weights.h5?download=1",
    "light":  "https://zenodo.org/records/13764771/files/light.weights.h5?download=1",
    "paired": "https://zenodo.org/records/13764771/files/paired.weights.h5?download=1",
}


# default params
PARAMS = [['CONV', 40, 10, 1],
          ['DROP', 0.2],
          ['POOL', 2, 1],
          ['FLAT'],
          ['DENSE', 300]]
ENCODING_DIM = 10   # kidera
SEQ_LEN = 200       # kasearch positions
PAD_LEN = 10        # padding between H and L chains


def create_cnn(units_per_layer, input_shape,
               activation, regularizer, out_dim=1):
    """
    Generate the CNN layers with a Keras wrapper.
    Code adapted from https://github.com/dahjan/DMS_opt

    Parameters
    ---
    units_per_layer: architecture features in list format, i.e.:
        Filter information: [CONV, # filters, kernel size, stride]
        Max Pool information: [POOL, pool size, stride]
        Dropout information: [DROP, dropout rate]
        Flatten: [FLAT]
        Dense layer: [DENSE, number nodes]

    input_shape: a tuple defining the input shape of the data

    activation: Activation function, i.e. ReLU, softmax

    regularizer: Kernel and bias regularizer in convulational and dense
        layers, i.e., regularizers.l1(0.01)
    """

    # Initialize the CNN
    model = keras.Sequential()

    # Input layer
    model.add(keras.layers.InputLayer(input_shape))

    # Build network
    for i, units in enumerate(units_per_layer):
        if units[0] == 'CONV':
            model.add(keras.layers.Conv1D(filters=units[1],
                                          kernel_size=units[2],
                                          strides=units[3],
                                          activation=activation,
                                          kernel_regularizer=regularizer,
                                          bias_regularizer=regularizer,
                                          padding='same'))
        elif units[0] == 'POOL':
            model.add(keras.layers.MaxPool1D(pool_size=units[1],
                                             strides=units[2]))
        elif units[0] == 'DENSE':
            model.add(keras.layers.Dense(units=units[1],
                                         activation=activation,
                                         kernel_regularizer=regularizer,
                                         bias_regularizer=regularizer))
        elif units[0] == 'DROP':
            model.add(keras.layers.Dropout(rate=units[1]))
        elif units[0] == 'FLAT':
            model.add(keras.layers.Flatten())
        else:
            raise NotImplementedError('Layer type not implemented')

    # Output layer
    # Activation function: Sigmoid if binary classification, softmax for multiclass
    activation = 'sigmoid' if out_dim == 1 else 'softmax'
    model.add(keras.layers.Dense(out_dim, activation=activation))

    return model


def load_cnn(weights, cnn_type, params=PARAMS):
    '''
    We save the checkpoint weights so need to load the relevant params too
    If retrained and full weights saved, use tf.keras.models.load_model(weights)

    :param weights: str, path to weights file
    :param cnn_type: str, type of CNN model heavy | light | paired
    :return: keras model
    '''
    if cnn_type == "heavy":
        seq_len, out_dim = SEQ_LEN, len(HEAVY_V_GENE_CLASSES)
        # check if weights file exists and download if not
        if not os.path.exists(weights):
            print(f"Downloading {cnn_type} model weights from {model_urls[cnn_type]}... ", end="")
            r = requests.get(model_urls[cnn_type], allow_redirects=True)
            open(weights, 'wb').write(r.content)
            print("done")
    elif cnn_type == "light":
        seq_len, out_dim = SEQ_LEN, len(LIGHT_V_GENE_CLASSES)
        if not os.path.exists(weights):
            print(f"Downloading {cnn_type} model weights from {model_urls[cnn_type]}... ", end="")
            r = requests.get(model_urls[cnn_type], allow_redirects=True)
            open(weights, 'wb').write(r.content)
            print("done")
    elif cnn_type == "paired":
        seq_len, out_dim = SEQ_LEN*2+PAD_LEN, len(PAIRED_CLASSES)
        if not os.path.exists(weights):
            print(f"Downloading {cnn_type} model weights from {model_urls[cnn_type]}... ", end="")
            r = requests.get(model_urls[cnn_type], allow_redirects=True)
            open(weights, 'wb').write(r.content)
            print("done")
    else:
        raise ValueError("cnn_type must be heavy | light | paired")
    
    CNN = create_cnn(params, (seq_len, ENCODING_DIM), 'relu', None, out_dim=out_dim)
    CNN.load_weights(weights)
    return CNN
