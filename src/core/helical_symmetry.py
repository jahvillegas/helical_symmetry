import numpy as np

def compute_helical_parameters(n1, n2, twist, rise):
    """
    Compute coordinates for helical assemblies based on the augmented 1-D helical system.
    Args:
        n1 (int): Number of helices.
        n2 (int): Translation part of the screw operation.
        twist (float): Twist angle in degrees.
        rise (float): Rise distance in Ångström.
    Returns:
        np.array: Array of transformed coordinates.
    """
    twist_rad = np.deg2rad(twist)
    transformation_matrix = np.array([
        [np.cos(twist_rad), -np.sin(twist_rad), 0],
        [np.sin(twist_rad),  np.cos(twist_rad), 0],
        [0, 0, 1]
    ])
    translation_vector = np.array([0, 0, rise])
    return transformation_matrix, translation_vector

