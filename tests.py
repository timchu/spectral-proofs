from spectral import *

# Tests. Can move to a different file.
# Would like some magic code that helps me draw graphs.
G = [[0, 1, 0], [1, 0, 0], [0, 0, 0]]
print(connected(G) == False)  # False
G = [[0, 1, 0], [1, 0, 1], [0, 1, 0]]
print(connected(G) == True)  # True
G = [[0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]]
print(connected(G) == False)  # False
G = [
    [0, 3, 0, 0, 0],
    [3, 0, 1, 0, 0],
    [0, 1, 0, 0.5, 0],
    [0, 0, 0.5, 0, 1],
    [0, 0, 0, 1, 0],
]
print(connected(G) == True)  # True
v = [5, 4, 3, 2, 1]
print(thres(G, v) == {3, 4})  # works as intended.
G = [
    [0, 3, 1, 0, 0, 0],
    [3, 0, 1, 0, 0, 0],
    [1, 1, 0, 1, 0, 0],
    [0, 0, 1, 0, 0, 0],
    [0, 0, 0, 0, 0, 1],
    [0, 0, 0, 0, 1, 0],
]
v = [6, 5, 4, 3, 2, 1]
print(connected(G) == False)
print(thres(G, v) == {4, 5})  # works as intended.
v = [1, 3, 2, 0, 5, 4]
print(thres(G, v) == {0, 1, 2, 3})  # works as intended.
G = [
    [0, 3, 1, 0, 0, 0],
    [3, 0, 1, 0, 0, 0],
    [1, 1, 0, 1, 0, 0],
    [0, 0, 1, 0, 0.1, 0],
    [0, 0, 0, 0.1, 0, 1],
    [0, 0, 0, 0, 1, 0],
]
v = [6, 5, 4, 3, 2, 1]
print(connected(G) == True)
print(thres(G, v) == {4, 5})  # works as intended.
v = [1, 3, 2, 0, 5, 4]
print(thres(G, v) == {0, 1, 2, 3})  # works as intended.
