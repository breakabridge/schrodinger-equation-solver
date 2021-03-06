import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

file = open("schrodinger_sim.dat")
data = file.readlines()
file.close()
line = data[0].split(",")

N = int(line[0])
frames = int(line[1])

dens = np.zeros((N, N), dtype = 'float')
line = data[1].split(",")
for i in range(N):
        for j in range(N):
            dens[i][j] = float(line[i + N * j])
fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on = True)
ax.set_xticks([])
ax.set_yticks([])
im = plt.imshow(dens, cmap = 'inferno')

def init():
    return []

def animate(t, im, dens):
    print(t)
    line = data[t+1].split(",")
    for i in range(N):
        for j in range(N):
            dens[i][j] = float(line[i + N * j])
    im.set_data(dens)
    return im,

if __name__ == "__main__":

    ani = animation.FuncAnimation(fig, animate, frames = frames, fargs = (im, dens), init_func = init, interval=500, blit=True)

    ani.save('video.mp4', fps = 24)
