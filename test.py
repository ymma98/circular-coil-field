import numpy as np
import matplotlib.pyplot as plt
from calcCurLoop import CoilArray


def calc_B():
    # generate coil array
    coil_num = 10
    coil_ra = 0.4 # radius of coil
    start_z = -2.
    end_z = 2.
    coil_current = 31.2e3
    coil_ra_list = np.ones(coil_num) * coil_ra
    coil_z_list = np.linspace(start_z, end_z, coil_num)
    coil_cur_list = np.ones(coil_num) * coil_current
    coil_arr = CoilArray(coil_ra_list,\
                         coil_z_list,\
                         coil_cur_list)
    # define target area
    rmin = 0.02
    rmax = 0.3
    zmin = -2.
    zmax = 2.
    mr = int((rmax-rmin)/0.01)
    mz = int((zmax-zmin)/0.05)
    r = np.linspace(rmin, rmax, mr)
    z = np.linspace(zmin, zmax, mz)
    rr, zz = np.meshgrid(r,z)
    # calculate br2d and bz2d in target area
    br2d = np.zeros((rr.shape[0], rr.shape[1]))
    bz2d = np.zeros((rr.shape[0], rr.shape[1]))
    for i in range(rr.shape[0]):
        for j in range(rr.shape[1]):
            br2d[i,j], bz2d[i,j] =\
                coil_arr.get_B(rr[i,j], zz[i,j])
    # save results
    np.savez(
           "./coil_mag.npz",
           coil_ra_list = coil_ra_list,
           coil_z_list = coil_z_list,
           coil_cur_list = coil_cur_list,
           r2d = rr,
           z2d = zz,
           br2d = br2d,
           bz2d = bz2d
       )


def plot_B(rr, zz, br2d, bz2d):
    # contour of br2d
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    cs1 = ax1.contourf(rr, zz, br2d, cmap='rainbow')
    fig1.colorbar(cs1)

    # contour of bz2d
    fig2 = plt.figure(2)
    ax2 = fig2.add_subplot(111)
    cs2 = ax2.contourf(rr, zz, bz2d, cmap="rainbow")
    fig2.colorbar(cs2)

    # streamline plot of magnetic fields
    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)
    ax3.streamplot(r2d, z2d, br2d, bz2d, cmap="rainbow", color=bz2d)

    plt.show()




if __name__=="__main__":
    # calculate B
    calc_B()

    ## plot data
    data = np.load("./coil_mag.npz")
    r2d = data['r2d']
    z2d = data['z2d']
    br2d = data['br2d']
    bz2d = data['bz2d']
    plot_B(r2d, z2d, br2d, bz2d)





