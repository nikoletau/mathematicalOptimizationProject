import time
from model import Model
from params import Params

p = Params()

start_time = time.time()
m = Model()
m.set_constraints(p)
m.set_obj(p)
m.optimise(logs=False)
end_time = time.time()
print("Computation time", (end_time - start_time))




