"""
ymma, ymma98@qq.com
2022.06.08
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import ellipk, ellipe


"""
ref:
[1] Jackson, Thrid edition, eq. (5.38)
[2] Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop
    James Simpson et al
notice: In Jackson's book, eq.(5.37) is wrong. For elliptic function on the numerator,
        they should be E(k**2) and K(k**2), but not the E(k) and K(k)

the derivations of Br and Bz are in ./sympy_derive.py
"""


class Coil:
    def __init__(self, r, z, cur):
        """
        :param r: (float) radius of the coil
        :param z: (float) location of the coil in z dir
        :param cur: (float) magnitude of the current
            counterclockwise: +, insert to the page
            clockwise: -, outward to the page
        """
        self.r = r
        self.z = z
        self.cur = cur


    def get_B(self, pr, pz):
        """
        :param pr: (float) coordinate r of (r, z) of given point
        :param pz: (float) coordinate z of (r, z) of given point
        """
        r = np.sqrt((0. - pr)**2 + (self.z - pz)**2) # distance between center of coil and the given point
        theta = np.arcsin(pr/r)
        Br_jackson = self._Br(pr, pz)
        Bth_jackson = self._Bth(pr,pz)
        # print("*****"*3)
        # print(r)
        # print(Br_jackson)
        # print(Bth_jackson)
        # print("*****"*3)
        if (pz > self.z):  # the point is above the coil
            Br = Bth_jackson * np.cos(theta) + Br_jackson * np.sin(theta)
            Bz = Br_jackson * np.cos(theta) - Bth_jackson * np.sin(theta)
        else:
            Br = -Bth_jackson * np.cos(theta) + Br_jackson * np.sin(theta)
            Bz = -Br_jackson * np.cos(theta) - Bth_jackson * np.sin(theta)
        return np.array([Br, Bz])


    def _Br(self, pr, pz):
        """
        :param pr: (float) coordinate r of (r, z) of given point
        :param pz: (float) coordinate z of (r, z) of given point
        """
        I = self.cur
        mu = 4*np.pi*1.e-7
        a = self.r
        r = np.sqrt((0. - pr)**2 + (self.z - pz)**2) # distance between center of coil and the given point
        if (pz<self.z):  # pz is under the coil
            theta = np.pi - np.arcsin(pr/r)
        else:
            theta = np.arcsin(pr/r)
        res = I*a**2*mu*np.cos(theta)*ellipe(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2))/(np.pi*(a**2 - 2*a*r*np.sin(theta) + r**2)*np.sqrt(a**2 + 2*a*r*np.sin(theta) + r**2))
        return res

    def _Bth(self, pr, pz):
        """
        :param pr: (float) coordinate r of (r, z) of given point
        :param pz: (float) coordinate z of (r, z) of given point
        """
        I = self.cur
        mu = 4*np.pi*1.e-7
        a = self.r
        r = np.sqrt((0. - pr)**2 + (self.z - pz)**2) # distance between center of coil and the given point
        if (pz<self.z):  # pz is under the coil
            theta = np.pi - np.arcsin(pr/r)
        else:
            theta = np.arcsin(pr/r)
        res =  I*mu*(a**2*np.cos(2*theta)*ellipe(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2)) - a**2*ellipk(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2)) + 2*a*r*np.sin(theta)*ellipk(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2)) + r**2*ellipe(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2)) - r**2*ellipk(4*a*r*np.sin(theta)/(a**2 + 2*a*r*np.sin(theta) + r**2)))/(2*np.pi*np.sqrt(a**2 + 2*a*r*np.sin(theta) + r**2)*(a**2*np.sin(theta) + a*r*np.cos(2*theta) - a*r + r**2*np.sin(theta)))
        return res



class CoilArray():
    def __init__(self, r_arr, z_arr, cur_arr):
        """
        :param r_arr: (np.ndarray) array of radii of coils
        :param z_arr: (np.ndarray) array of locations in z dir of coils
        :param cur_arr: (np.ndarray) array of current of coils
            counterclockwise: +, insert to the page
            clockwise: -, outward to the page
        """
        self.num = r_arr.shape[0]  # number of coils
        self.r_arr = r_arr
        self.z_arr = z_arr
        self.cur_arr = cur_arr
        self.coil_arr = []
        for i in range(self.num):
            self.coil_arr.append( \
                Coil(self.r_arr[i],\
                     self.z_arr[i],\
                     self.cur_arr[i]))

    def get_B(self, pr, pz):
        """
        :param pr: (float) coordinate r of (r, z) of given point
        :param pz: (float) coordinate z of (r, z) of given point
        """
        br = 0.
        bz = 0.
        for i in range(self.num):
            tmpbr,tmpbz = self.coil_arr[i].get_B(pr, pz)
            br += tmpbr
            bz += tmpbz
        return (br, bz)


