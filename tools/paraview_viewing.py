import afidtools as afid
import numpy as np
folder = "/scratch/seismo/dave/melting_Ross/0.3"
vars = ["temp","phi"]
for var in vars:
    if var == "phi":
        afid.interpolate_field_to_uniform(folder, var, scale=2)
        print("Interpolated "+var+" to uniform")
        afid.generate_uniform_xmf(folder, var, scale=2)
        print("Created xmf for "+var)
    else:
        afid.interpolate_field_to_uniform(folder, var, scale=1)
        print("Interpolated "+var+" to uniform")
        afid.generate_uniform_xmf(folder, var, scale=1)
        print("Created xmf for "+var)


    