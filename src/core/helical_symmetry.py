import numpy as np

def compute_helical_parameters(twist, rise):
    twist_rad = np.deg2rad(twist)
    #twist_rad = twist
    rotation_matrix = np.array([
        [np.cos(twist_rad), -np.sin(twist_rad), 0],
        [np.sin(twist_rad),  np.cos(twist_rad), 0],
        [0,                  0,                 1]
    ])
    translation_vector = np.array([0, 0, rise])
    return rotation_matrix, translation_vector
