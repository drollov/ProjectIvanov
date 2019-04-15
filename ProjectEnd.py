import math
import time
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
#%%
def len_of_time():
    return print("--- %s seconds ---" % (time.time() - start_time)) 

def graph_bandwidth(users, bandwidth, S):
    color = colors.pop()
    R_max = "%.3f" % max(bandwidth)
#    plt.plot(users[1:], bandwidth[1:], color = color, linewidth = 2) 
    plt.plot(users, bandwidth, color = color, label = '(S_theor, R_max) = ' + str((S, R_max)), linewidth = 2) 
    
#    handles.append(mpatches.Patch(color = color, label = 'S = ' + str(S)))
    plt.legend(loc = 'lower left')
#    plt.legend(handles=handles)
    plt.xlabel('U', rotation = 0)
    plt.ylabel('R_sigma', rotation = 0)
#    plt.ylabel('C * L', rotation = 0)
#    plt.title("""Зависимость между T и числом активных пользователей U""")
    plt.grid(True)
    plt.show()
    
    print('Max of bandwidth = ' + str(max(bandwidth)))
    
def combination(length, w):
    result = 1
    denomination = 1
    while w > 0:
        result *= (length - (w - 1)) / denomination
        denomination += 1
        w -= 1
    return result
#%%
class Analysis:
    def __init__(self, U, S, L):
        self.U = U #Число активных пользователей
        self.S = S #Число независимых каналов
        self.L = L #Кол-во слотов в независимом канале
    
    def P_prob_k(self, k):
        return combination(self.U - 1, k) * (1/self.S)**k * (1 - 1/self.S)**(self.U - k - 1)
    
    def P_in(self, w):
        return combination(self.L, w) * 2**(-self.L)
    
    def p(self, k):
        return 1 - 2**(-k)
    
    def P(self, w_check, w, k):
        if w_check < w:
            return 0
        else:
            p_k = self.p(k)
            return combination(self.L - w, w_check - w) * \
                   p_k**(w_check - w) *(1 - p_k)**(self.L - w_check)
                       
    def P_out(self, w_check, k):
        Sum = 0
        for i in range(w_check + 1):
            Sum += self.P(w_check, i, k) * self.P_in(i)
        return Sum
    
    def min_bandwidth_collision(self, k):
        outer_sum = 0
        minus_sum = 0
        
        for i in range(self.L + 1):
            P_out = self.P_out(i, k)
            if P_out > 0:
                minus_sum += P_out * math.log2(P_out)
                
            inner_sum = 0
            for j in range(i, self.L + 1):
                P = self.P(j, i, k)
                if P > 0:
                    inner_sum += P * math.log2(P)
            outer_sum += self.P_in(i) * inner_sum          
         
        return outer_sum - minus_sum
    
    def band_coll(self):
        Sum = 0
        print(self.U)
        for i in range(self.U):
            if i == 0:
                Sum += self.P_prob_k(i)
            else:
                Sum += self.P_prob_k(i) * self.min_bandwidth_collision(i)
        
        return Sum
#%%
class VG_fix_delta:
    def __init__(self, L, Q, U, delta):
        self.L = L
        self.Q = Q
        self.U = U
        self.q = 2 ** self.L
        self.S = self.Q / self.L
        self.delta = delta
        
    def P_s(self, qty):
        return 1 - (1 - 1 / self.S)**(qty - 1)
    
    def VG(self):   
        asymp = -self.delta * math.log(self.delta, self.q) - (1 - self.delta) * math.log(1 - self.delta, self.q) + \
        self.delta * math.log(self.q - 1, self.q)
        
        return 1 - asymp
    
    def binom_distrib(self, qty, n):
        Sum = 0
        for i in np.arange(self.delta * n):
            Sum += combination(n, i) * pow(self.P_s(qty), i) * pow(1 - self.P_s(qty), n - i) 
        
        return Sum
           
    def graph_fix_delta(self, n):
#        color = colors.pop()
        
        mas = [i * self.L * self.VG() * self.binom_distrib(i, n) / self.Q for i in range(self.U)]
        r_max = "%.3f" % max(mas)
        extrema.append((max(mas), n))
        plt.plot(list(range(self.U)), mas, colors.pop(), label = '(S_VG, R_max) = ' + str((self.S, r_max)))
        plt.legend(loc = 'lower right')
#        handles.append(mpatches.Patch(color = color, label = 'n = ' + str(n)))

#        plt.legend(handles=handles)
#        
#        plt.xlabel('U', rotation = 0)
#        plt.ylabel('R_sigma', rotation = 0)
        plt.grid(True)
        plt.show() 
#%%            
if __name__ == "__main__":
    n = 170
    delta = 0.6
    qty_users = int(input('Кол-во пользователей = ')) #желательно для 500
    start_time = time.time()
    colors = ['g', 'r', 'b', 'm', 'y', 'c', '#7B2424', '#3C7B24']
    handles = []
    extrema = []
    L = 128

    for Q in [512, 1024, 2048]:
    
        S = Q / L 
        print('S = ' + str(S))
                        
        users = list(range(qty_users))
        
        #        bandwidth = list(map(lambda U: Analysis(U, S, L).band_coll() * L, users))      
        #        graph_bandwidth(users, bandwidth, S)
        
        R_sigma = list(map(lambda U: Analysis(U, S, L).band_coll() * U * L / Q, users)) #Спектральная эффективность (она показывает относительный объем инфы, которые могут пердеать U пользователей)       
        graph_bandwidth(users, R_sigma, S)
            
        print(n)
        VG_fix_delta(L, Q, len(users), delta).graph_fix_delta(n)            
    
    print(max(extrema))
    len_of_time()
        
#     Кодовое слово, в случае если возникли коллизии 
#    n - число слотов в кодовом слове, k - число информац словтов, d - расстояние между двумя словами 
#  код задан над конечным полем!!    
#       
        
        
        
        
    
