import pdb

import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.python.keras import backend as K
from tensorflow.python.keras.utils import tf_utils
from tensorflow.python.ops import array_ops


class Encoder_T(keras.layers.Layer):
    """
    Encoder for transcriptomic data
    
    Args:
        dropout_rate: dropout probability if training=True
        latent_dim: dimensionality of representation
        intermediate_dim: number of units in hidden layers
    """

    def __init__(self,
                 dropout_rate=0.5,
                 latent_dim=3,
                 intermediate_dim=50,
                 name='Encoder_T',
                 dtype=tf.float32,
                 **kwargs):

        super(Encoder_T, self).__init__(name=name, **kwargs)
        self.drp = keras.layers.Dropout(rate=dropout_rate)
        self.fc0 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc0')
        self.fc1 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc1')
        self.fc2 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc2')
        self.fc3 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc3')
        self.fc4 = keras.layers.Dense(latent_dim, activation='linear', name=name+'fc4')
        self.bn = keras.layers.BatchNormalization(scale=False, center=False, epsilon=1e-10, momentum=0.05, name=name+'BN')
        return

    def call(self, inputs, training=True):
        x = self.drp(inputs, training=training)
        x = self.fc0(x, training=training)
        x = self.fc1(x, training=training)
        x = self.fc2(x, training=training)
        x = self.fc3(x, training=training)
        x = self.fc4(x, training=training)
        z = self.bn(x, training=training)
        return z


class Decoder_T(keras.layers.Layer):
    """
    Decoder for transcriptomic data

    Args:
        output_dim: number of outputs
        intermediate_dim: number of units in hidden layers
    """

    def __init__(self,
                 output_dim,
                 intermediate_dim=50,
                 name='Decoder_T',
                 dtype=tf.float32,
                 **kwargs):
        
        super(Decoder_T, self).__init__(name=name, **kwargs)
        self.fc0 = keras.layers.Dense(intermediate_dim, activation='relu', name='fc0')
        self.fc1 = keras.layers.Dense(intermediate_dim, activation='relu', name='fc1')
        self.fc2 = keras.layers.Dense(intermediate_dim, activation='relu', name='fc2')
        self.fc3 = keras.layers.Dense(intermediate_dim, activation='relu', name='fc3')
        # transcriptomics data should be larger than 0?
        self.Xout = keras.layers.Dense(output_dim, activation='relu', name='Xout')
        return

    def call(self, inputs, training=True):
        x = self.fc0(inputs, training=training)
        x = self.fc1(x, training=training)
        x = self.fc2(x, training=training)
        x = self.fc3(x, training=training)
        x = self.Xout(x)
        return x


class Encoder_E(keras.layers.Layer):
    """
    Decoder for electrophysiology data
    
    Args:
        gaussian_noise_sd: std of gaussian noise injection if training=True
        dropout_rate: dropout probability if training=True
        latent_dim: representation dimenionality
        intermediate_dim: number of units in hidden layers
    """

    def __init__(self,
                 gaussian_noise_sd=0.05,
                 dropout_rate=0.1,
                 latent_dim=3,
                 intermediate_dim=40,
                 name='Encoder_E',
                 dtype=tf.float32,
                 **kwargs):
        
        super(Encoder_E, self).__init__(name=name, **kwargs)
        self.gnoise = WeightedGaussianNoise(stddev=gaussian_noise_sd)
        self.drp = keras.layers.Dropout(rate=dropout_rate)
        self.fc0 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc0')
        self.fc1 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc1')
        self.fc2 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc2')
        self.fc3 = keras.layers.Dense(intermediate_dim, activation='relu', name=name+'fc3')
        self.fc4 = keras.layers.Dense(latent_dim, activation='linear', name=name+'fc4')
        self.bn = keras.layers.BatchNormalization(scale=False, center=False, epsilon=1e-10, momentum=0.05, name=name+'BN')
        return

    def call(self, inputs, training=True):
        x = self.gnoise(inputs, training=training)
        x = self.drp(x, training=training)
        x = self.fc0(x, training=training)
        x = self.fc1(x, training=training)
        x = self.fc2(x, training=training)
        x = self.fc3(x, training=training)
        x = self.fc4(x, training=training)
        z = self.bn(x, training=training)
        return z

class Decoder_E(keras.layers.Layer):
    """
    Initializes the Encoder for electrophysiology data.

    Args:
        output_dim: Should be same as input dim if using as an autoencoder
        intermediate_dim: Number of units in hidden keras.layers
        training: boolean value to indicate model operation mode
    """

    def __init__(self,
                 output_dim,
                 intermediate_dim=40,
                 name='Decoder_E',
                 dtype=tf.float32,
                 **kwargs):
   
        super(Decoder_E, self).__init__(name=name, **kwargs)
        self.fc0 = keras.layers.Dense(intermediate_dim, activation='relu',name=name+'fc0')
        self.fc1 = keras.layers.Dense(intermediate_dim, activation='relu',name=name+'fc1')
        self.fc2 = keras.layers.Dense(intermediate_dim, activation='relu',name=name+'fc2')
        self.fc3 = keras.layers.Dense(intermediate_dim, activation='relu',name=name+'fc3')
        self.drp = keras.layers.Dropout(rate=0.1) 
        self.Xout = keras.layers.Dense(output_dim, activation='linear',name=name+'Xout')
        return

    def call(self, inputs, training=True):
        x = self.fc0(inputs, training=training)
        x = self.fc1(x, training=training)
        x = self.fc2(x, training=training)
        x = self.fc3(x, training=training)
        x = self.drp(x, training=training)
        x = self.Xout(x, training=training)
        return x


class Model_TE(tf.keras.Model):
    """
    Coupled autoencoder

    Args:
        T_dim: n(genes)
        E_dim: n(features)
        T_intermediate_dim: units in hidden layers of T autoencoder
        E_intermediate_dim: units in hidden layers of E autoencoder
        T_dropout: dropout probability for 
        E_gnoise_sd: gaussian noise std for E data
        E_dropout: dropout for E data
        latent_dim: dim for representations
        name: TE

    call Args:
        train_T: augment T data
        train_E: augment E data
    
    """

    def __init__(self,
               T_dim,
               E_dim,
               T_intermediate_dim=50,
               E_intermediate_dim=40,
               T_dropout=0.5,
               E_gnoise_sd=0.05,
               E_dropout=0.1,
               latent_dim=3,
               name='TE',
               **kwargs):
  
        super(Model_TE, self).__init__(name=name, **kwargs)
        self.encoder_T = Encoder_T(dropout_rate=T_dropout,
                                   latent_dim=latent_dim,
                                   intermediate_dim=T_intermediate_dim,
                                   name='Encoder_T')

        self.encoder_E = Encoder_E(gaussian_noise_sd=E_gnoise_sd,
                                   dropout_rate=E_dropout,
                                   latent_dim=latent_dim,
                                   intermediate_dim=E_intermediate_dim,
                                   name='Encoder_E')

        self.decoder_T = Decoder_T(output_dim=T_dim,
                                   intermediate_dim=T_intermediate_dim,
                                   name='Decoder_T')

        self.decoder_E = Decoder_E(output_dim=E_dim,
                                   intermediate_dim=E_intermediate_dim,
                                   name='Decoder_E')

    def call(self, inputs, train_T=True, train_E=True):
        #T arm
        zT = self.encoder_T(inputs[0],training=train_T)
        zE = self.encoder_E(inputs[1],training=train_E)

        XrT = self.decoder_T(zT,training=train_T)
        XrE = self.decoder_E(zE,training=train_E)
        
        XrT_cm = self.decoder_T(zE,training=train_E)
        XrE_cm = self.decoder_E(zT,training=train_T)
        return zT,zE,XrT,XrE,XrT_cm,XrE_cm


class WeightedGaussianNoise(keras.layers.Layer):
    """Custom additive zero-centered Gaussian noise. Std is weighted.

    Args:
        stddev: Can be a scalar or vector
    call args:
        inputs: Input tensor (of any rank).
        training: Python boolean indicating whether the layer should behave in
        training mode (adding noise) or in inference mode (doing nothing).
    Input shape:
        Arbitrary. Use the keyword argument `input_shape`
        (tuple of integers, does not include the samples axis)
        when using this layer as the first layer in a model.
    Output shape:
        Same shape as input.
    """

    def __init__(self, stddev, **kwargs):
        super(WeightedGaussianNoise, self).__init__(**kwargs)
        self.stddev = stddev
        return

    def call(self, inputs, training=None):
        def noised():
            return inputs + tf.random.normal(array_ops.shape(inputs),
                                             mean=0.0, stddev=self.stddev,
                                             dtype=inputs.dtype, seed=None)

        return K.in_train_phase(noised, inputs, training=training)

    def get_config(self):
        config = {'stddev': self.stddev}
        base_config = super(WeightedGaussianNoise, self).get_config()
        return dict(list(base_config.items()) + list(config.items()))

    @tf_utils.shape_type_conversion
    def compute_output_shape(self, input_shape):
        return input_shape


class Model_TE_aug_decoders(tf.keras.Model):
    """Coupled autoencoder model

    Args:
        T_dim: Number of genes in T data
        E_dim: Number of features in E data
        T_intermediate_dim: hidden layer dims for T model
        E_intermediate_dim: hidden layer dims for E model
        T_dropout: dropout for T data
        E_gnoise_sd: gaussian noise std for E data
        E_dropout: dropout for E data
        latent_dim: dim for representations
        name: TE
    """
    def __init__(self,
               T_dim,
               E_dim,
               T_intermediate_dim=50,
               E_intermediate_dim=40,
               alpha_T=1.0,
               alpha_E=1.0,
               lambda_TE=1.0,
               T_dropout=0.5,
               E_gauss_noise_wt = 1.0,
               E_gnoise_sd=0.05,
               E_dropout=0.1,
               latent_dim=3,
               name='TE',
               **kwargs):

        super(Model_TE_aug_decoders, self).__init__(name=name, **kwargs)
        self.T_dim = T_dim
        self.E_dim = E_dim

        self.alpha_T = tf.constant(alpha_T,dtype=tf.float32)
        self.alpha_E = tf.constant(alpha_E,dtype=tf.float32)
        self.lambda_TE = tf.constant(lambda_TE,dtype=tf.float32)

        E_gnoise_sd_weighted = E_gauss_noise_wt*E_gnoise_sd
        self.encoder_T = Encoder_T(dropout_rate=T_dropout,latent_dim=latent_dim, intermediate_dim=T_intermediate_dim, name='Encoder_T')
        self.encoder_E = Encoder_E(gaussian_noise_sd=E_gnoise_sd_weighted, dropout_rate=E_dropout, latent_dim=latent_dim, intermediate_dim=E_intermediate_dim, name='Encoder_E')
        
        self.decoder_T = Decoder_T(output_dim=T_dim, intermediate_dim=T_intermediate_dim, name='Decoder_T')
        self.decoder_E = Decoder_E(output_dim=E_dim, intermediate_dim=E_intermediate_dim, name='Decoder_E')
        return


    def call(self, inputs, train_T=True, train_E=True, augment_decoders=True):
        """
        Args:
            train_T: training/inference mode for T autoencoder
            train_E: training/inference mode for E autoencoder
            augment_decoders: augment decoder with cross modal representation if True
        """
        #T arm forward pass
        XT = inputs[0]
        zT = self.encoder_T(XT,training=train_T)
        XrT = self.decoder_T(zT,training=train_T)
        
        #E arm forward pass
        XE = tf.where(tf.math.is_nan(inputs[1]),x=0.0,y=inputs[1]) #Mask nans
        maskE = tf.where(tf.math.is_nan(inputs[1]),x=0.0,y=1.0)    #Get mask to ignore error contribution
        zE = self.encoder_E(XE,training=train_E)
        XrE = self.decoder_E(zE,training=train_E)

        #Loss calculations
        mse_loss_T = tf.reduce_mean(tf.math.squared_difference(XT, XrT))
        mse_loss_E = tf.reduce_mean(tf.multiply(tf.math.squared_difference(XE, XrE),maskE))
        cpl_loss_TE = min_var_loss(zT, zE)

        #Append to keras model losses for gradient calculations
        self.add_loss(self.alpha_T*mse_loss_T)
        self.add_loss(self.alpha_E*mse_loss_E)
        self.add_loss(self.lambda_TE*cpl_loss_TE)

        #Cross modal reconstructions - treat zE and zT as constants for this purpose
        if augment_decoders:
            XrT_aug = self.decoder_T(tf.stop_gradient(zE),training=train_T)
            XrE_aug = self.decoder_E(tf.stop_gradient(zT),training=train_E)
            mse_loss_T_aug = tf.reduce_mean(tf.math.squared_difference(XT, XrT_aug))
            mse_loss_E_aug = tf.reduce_mean(tf.multiply(tf.math.squared_difference(XE, XrE_aug),maskE))
            self.add_loss(self.alpha_T*mse_loss_T_aug)
            self.add_loss(self.alpha_E*mse_loss_E_aug)
        
        #For logging only
        self.mse_loss_T = mse_loss_T
        self.mse_loss_E = mse_loss_E
        self.mse_loss_TE = tf.reduce_mean(tf.math.squared_difference(zT, zE))
        self.mse_loss_T_aug = mse_loss_T_aug
        self.mse_loss_E_aug = mse_loss_E_aug
        
        XrT_cm = self.decoder_T(zE,training=train_E)
        XrE_cm = self.decoder_E(zT,training=train_T)
        return zT,zE,XrT,XrE,XrT_cm,XrE_cm
        

def min_var_loss(zi, zj, Wij=None):
    """
    SVD is calculated over entire batch. MSE is calculated over only paired entries within batch
    
    Args:
        zi: i-th representation (batch_size x latent_dim)
        zj: j-th representation (batch_size x latent_dim)
        Wij: indicator vector (batch_size x latent_dim) (1 if samples are paired, 0 otherwise)
    """
    batch_size = tf.shape(zi)[0]
    if Wij is None:
        Wij_ = tf.ones([batch_size, ])
    else:
        Wij_ = tf.reshape(Wij, [batch_size, ])

    #Masking gets rid of unpaired entries
    zi_paired = tf.boolean_mask(zi, tf.math.greater(Wij_, 1e-2))
    zj_paired = tf.boolean_mask(zj, tf.math.greater(Wij_, 1e-2))
    Wij_paired = tf.boolean_mask(Wij_, tf.math.greater(Wij_, 1e-2))

    #SVD calculated over all entries in the batch
    vars_j_ = tf.square(tf.reduce_min(tf.linalg.svd(zj - tf.reduce_mean(zj, axis=0), compute_uv=False)))/tf.cast(batch_size - 1, tf.float32)
    vars_j  = tf.where(tf.math.is_nan(vars_j_), tf.cast(1e-2,dtype=tf.float32), vars_j_)

    vars_i_ = tf.square(tf.reduce_min(tf.linalg.svd(zi - tf.reduce_mean(zi, axis=0), compute_uv=False)))/tf.cast(batch_size - 1, tf.float32)
    vars_i  = tf.where(tf.math.is_nan(vars_i_), tf.cast(1e-2,dtype=tf.float32), vars_i_)

    #Wij_paired is the weight of matched pairs
    sqdist_paired = tf.multiply(tf.reduce_sum(tf.math.squared_difference(zi_paired, zj_paired),axis=1),Wij_paired)
    
    mean_sqdist = tf.reduce_sum(sqdist_paired,axis=None)/tf.reduce_sum(Wij_paired,axis=None)
    loss_ij = mean_sqdist/tf.maximum(tf.reduce_min([vars_i,vars_j], axis=None),tf.cast(1e-2,dtype=tf.float32))
    return loss_ij


def custom_build(model,input_shape=None):
    """
    Initialize the network using this if loading saved weights. 
    
    Args: 
        input_shape: (T_dim, E_dim) tuple with input shapes for the two arms 
    """
    x = tf.constant(np.random.rand(1,input_shape[0]),dtype=tf.float32)
    y = tf.constant(np.random.rand(1,input_shape[1]),dtype=tf.float32)
    _,_,_,_,_,_ = model((x,y),train_T=False,train_E=False)
    return model