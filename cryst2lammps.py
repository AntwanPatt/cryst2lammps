import numpy as np

class cryst2lammps:

    def __init__(self, crystfile=None):
        self.infile = crystfile

    def set_base(self, edges=None, angles=90., angle_in_degrees=True):
        ## form base matrix A
        ## converting crystallographic vectors in cartesian coordinates
        ## using a formula from: (got at Chemistry StackExchange)

        try:
            a, b, c = edges
        except:
            a = b = c = edges

        try:
            alpha, beta, gamma = angles
        except:
            alpha = beta = gamma = angles

        if angle_in_degrees:
            alpha = np.deg2rad(alpha)
            beta = np.deg2rad(beta)
            gamma = np.deg2rad(gamma)

        n2 = (np.cos(alpha) - np.cos(gamma) * np.cos(beta)) / np.sin(gamma)
        factz = np.sqrt(np.sin(beta)**2. - n2**2)

        A = np.array([[a,                 0.,                0.        ],
                      [b * np.cos(gamma), b * np.sin(gamma), 0.        ],
                      [c * np.cos(beta),  c * n2,            c * factz ]])

        cosa = np.cos(alpha)
        sina = np.sin(alpha)
        cosb = np.cos(beta)
        sinb = np.sin(beta)
        cosg = np.cos(gamma)
        sing = np.sin(gamma)
        volume = 1.0 - cosa**2.0 - cosb**2.0 - cosg**2.0 + 2.0 * cosa * cosb * cosg
        volume = np.sqrt(volume)

        r = np.zeros((3, 3))
        r[0, 0] = a
        r[0, 1] = b * cosg
        r[0, 2] = c * cosb
        r[1, 1] = b * sing
        r[1, 2] = c * (cosa - cosb * cosg) / sing
        r[2, 2] = c * volume / sing

        return r


    #def read(self):
        ## read initial file for the specified format
        ## cartesian or direct coordinates in the file?

    #def replicate(self):
        ## multiply it
        ## condition: if bonds already set up can't replicate

    #def bonds(self):
        ## finding bonds
