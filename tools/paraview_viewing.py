import afidtools as afid
import numpy as np
folder = "/scratch/seismo/dave/melting_Ross/0.3"
vars = ["temp"]
for var in vars:
    afid.interpolate_field_to_uniform(folder, var, scale=4)
    print("Interpolated "+var+" to uniform")
    afid.generate_uniform_xmf(folder, var, scale=4)
    print("Created xmf for "+var)