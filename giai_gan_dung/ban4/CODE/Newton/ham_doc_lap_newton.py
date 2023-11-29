import sympy as sym
import scipy as sci
import numpy as np
from math import *
import sys


#===================================================================================
# Phần thuật toán chính
class newton_mat_inversion:
    def __init__(self, A, n, eps):    # Khởi tạo
        self.A = np.reshape(A, (n, n))
        self.n = n
        self.eps = eps
        self.nr_iterations = 0

    def norm(self, __A, __norm_type = 2):   # Chuẩn ma trận
        return np.linalg.norm(__A, __norm_type)

    def __refine_initial_approx(self):  # Tìm xấp xỉ đầu:
        # Gán các biến cơ bản
        E  = np.eye(self.n)
        A  = self.A
        X  = self.A

        t1 = self.norm(X, 1)
        t2 = self.norm(X, inf)
        X = (X / (t1 * t2)).T
        print("Ma trận X0:")
        print(X)
        return X

    def __pure_newton(self, X_0):  # Lặp Newton
        # Gán các biến cơ bản
        norm_X0 = self.norm(X_0) # Do X_0 không đổi nên ta đặt 1 biến làm chuẩn của X_0
        E       = np.eye(self.n)
        A       = self.A
        eps     = self.eps
        # Bước 3 của thuật toán
        q2k = q = self.norm(E - A @ X_0)
        X = X_0
        # Kiểm tra điều kiện hội tụ
        if(q >= 1):
            print("Xấp xỉ đầu không thỏa mãn nên không đưa ra được ma trận nghịch đảo.")
            return np.full((self.n, self.n), float("NaN"))
        # Lặp
        while(norm_X0 * q2k > self.eps * (1 - q)):
            self.nr_iterations += 1
            X = X @ (2 * E - A @ X)
            q2k = q2k ** 2
        # Đưa ra ma trận cuối cùng
        print(f"Phương pháp Newton kết thúc sau {self.nr_iterations} bước lặp", file=sys.stderr)
        return X

    def approx_newton(self):  # PP Newton với xấp xỉ đầu
        if(np.linalg.det(self.A) == 0):
            print("A không khả nghịch nên không đưa ra được ma trận nghịch đảo")
            return np.full((self.n, self.n), float("NaN"))
        return self.__pure_newton(self.__refine_initial_approx())


