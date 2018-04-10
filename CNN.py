import numpy as np
from keras.models import load_model

noise_hat_matrix = []
noise_wave_matrix = np.zeros((1000, 576))
noise_hat_matrix = np.loadtxt('noise_hat_t0.5_s0.txt')
model = load_model('model_0.5.h5')
for i in range(1000):
    print(i)
    noise_hat = noise_hat_matrix[i,:]
    noise_hat = noise_hat.reshape(1,576)
    noise_hat_reshape = np.expand_dims(noise_hat, axis=2)
    noise_wave = model.predict(noise_hat_reshape)
    noise_wave_reshape = noise_wave.reshape((1,576))
    noise_wave_matrix[i,:] = noise_wave_reshape 

np.savetxt("noise_wave_t0.5_s0.txt",noise_wave_matrix)