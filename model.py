import gurobipy as gp
from params import Params

class Model:
    def __init__(self):
        self.par = Params()
        self.m = gp.Model("salmone")
        self.transportation_cost = 0
        self.fuel_cost = 0
        self.residual_cost = 0
        self.inventory_cost = 0
        self.processing_cost = 0
        self.carbon_emission = 0

        # CONTINUOUS VARIABLES related to the processed amount & wastage amount
        # total processed amount of HOG 'b' at slaughterhouse 'j'
        self.tp_b = self.m.addVars([(j, t) for j in self.par.J for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for j in self.par.J:
            for t in self.par.T:
                self.tp_b[j, t].setAttr('VarName', f"tp_b{j}{t}")
        # total processed amount of fresh HOG 'c' at prim. proc. pl. 'p'
        self.tp_c = self.m.addVars([(p, t) for p in self.par.P for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for p in self.par.P:
            for t in self.par.T:
                self.tp_c[p, t].setAttr('VarName', f"tp_c{p}{t}")
        # total processed amount of fillet 'e' at sec. proc. pl. 'q'
        self.tp_e = self.m.addVars([(q, t) for q in self.par.Q for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for q in self.par.Q:
            for t in self.par.T:
                self.tp_e[q, t].setAttr('VarName', f"tp_e{q}{t}")
        # total processed amount of by-product 'f' at sec. proc. pl. 'q'
        self.tp_f = self.m.addVars([(q, t) for q in self.par.Q for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for q in self.par.Q:
            for t in self.par.T:
                self.tp_f[q, t].setAttr('VarName', f"tp_f{q}{t}")
        # amount of residual obtained after processing live salmon 'a' at slaughterhouse 'j'
        self.tw_a = self.m.addVars([(j, t) for j in self.par.J for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for j in self.par.J:
            for t in self.par.T:
                self.tw_a[j, t].setAttr('VarName', f"tw_a{j}{t}")
        # amount of residual obtained after processing HOG 'b' at prim. proc. pl. 'p'
        self.tw_b = self.m.addVars([(p, t) for p in self.par.P for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for p in self.par.P:
            for t in self.par.T:
                self.tw_b[p, t].setAttr('VarName', f"tw_b{p}{t}")
        # amount of residual obtained after processing fresh HOG 'c' at sec. proc. pl. 'q'
        self.tw_c = self.m.addVars([(q, t) for q in self.par.Q for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for q in self.par.Q:
            for t in self.par.T:
                self.tw_c[q, t].setAttr('VarName', f"tw_c{q}{t}")

        # CONTINUOUS VARIABLES related to the inventory level
        # inventory of HOG 'b' at slaughterhouse 'j'
        self.ip_b = self.m.addVars([(j, t) for j in self.par.J for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for j in self.par.J:
            for t in self.par.T:
                self.ip_b[j, t].setAttr('VarName', f"ip_b{j}{t}")
        # inventory of fresh HOG 'c' at prim. proc. pl. 'p'
        self.ip_cp = self.m.addVars([(p, t) for p in self.par.P for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for p in self.par.P:
            for t in self.par.T:
                self.ip_cp[p, t].setAttr('VarName', f"ip_cp{p}{t}")
        # inventory of fresh HOG 'c' at wholesaler 'u'
        self.ip_cu = self.m.addVars([(u, t) for u in self.par.U for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for u in self.par.U:
            for t in self.par.T:
                self.ip_cu[u, t].setAttr('VarName', f"ip_cu{u}{t}")
        # inventory of fillet 'e' at sec. proc. pl. 'q'
        self.ip_eq = self.m.addVars([(q, t) for q in self.par.Q for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for q in self.par.Q:
            for t in self.par.T:
                self.ip_eq[q, t].setAttr('VarName', f"ip_eq{q}{t}")
        # inventory of by-product 'f' at sec. proc. pl. 'q'
        self.ip_fq = self.m.addVars([(q, t) for q in self.par.Q for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for q in self.par.Q:
            for t in self.par.T:
                self.ip_fq[q, t].setAttr('VarName', f"ip_fq{q}{t}")
        # inventory of fillet 'e' at wholesaler 'u'
        self.ip_eu = self.m.addVars([(u, t) for u in self.par.U for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for u in self.par.U:
            for t in self.par.T:
                self.ip_eu[u, t].setAttr('VarName', f"ip_eu{u}{t}")
        # inventory of by-product 'f' at wholesaler 'u'
        self.ip_fu = self.m.addVars([(u, t) for u in self.par.U for t in self.par.T], vtype=gp.GRB.CONTINUOUS)
        for u in self.par.U:
            for t in self.par.T:
                self.ip_fu[u, t].setAttr('VarName', f"ip_fu{u}{t}")

       # INTEGER VARIABLES related to amount transported
        # total amount of live salmon 'a' transported from farm 'i' to slaughterhouse 'j'
        self.x_a = self.m.addVars([(i, j, t) for i in self.par.I for j in self.par.J for t in self.par.T], vtype=gp.GRB.INTEGER)
        for i in self.par.I:
                for j in self.par.J:
                    for t in self.par.T:
                        self.x_a[i, j, t].setAttr('VarName', f"x_a{i}{j}{t}")
        # total amount of HOG 'b' transported from slaughterhouse 'j' to prim. proc. pl. 'p'
        self.x_b = self.m.addVars([(j, p, t) for j in self.par.J for p in self.par.P for t in self.par.T], vtype=gp.GRB.INTEGER)
        for j in self.par.J:
            for p in self.par.P:
                for t in self.par.T:
                    self.x_b[j, p, t].setAttr('VarName', f"x_b{j}{p}{t}")
        # total amount of fresh HOG 'c' transported from prim. proc. pl. 'p' to sec. proc. pl. 'q'
        self.x_cpq = self.m.addVars([(p, q, t) for p in self.par.P for q in self.par.Q for t in self.par.T], vtype=gp.GRB.INTEGER)
        for p in self.par.P:
            for q in self.par.Q:
                for t in self.par.T:
                    self.x_cpq[p, q, t].setAttr('VarName', f"x_cpq{p}{q}{t}")
        # total amount of fresh HOG 'c' transported from prim. proc. pl. 'p' to wholesaler 'u'
        self.x_cpu = self.m.addVars([(p, u, t) for p in self.par.P for u in self.par.U for t in self.par.T], vtype=gp.GRB.INTEGER)
        for p in self.par.P:
            for u in self.par.U:
                for t in self.par.T:
                    self.x_cpu[p, u, t].setAttr('VarName', f"x_cpu{p}{u}{t}")
        # total amount of fillet 'e' transported from sec. proc. pl. 'q' to wholesaler 'u'
        self.x_equ = self.m.addVars([(q, u, t) for q in self.par.Q for u in self.par.U for t in self.par.T], vtype=gp.GRB.INTEGER)
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T:
                    self.x_equ[q, u, t].setAttr('VarName', f"x_equ{q}{u}{t}")
        # total amount of by-product 'f' transported from sec. proc. pl. 'q' to wholesaler 'u'
        self.x_fqu = self.m.addVars([(q, u, t) for q in self.par.Q for u in self.par.U for t in self.par.T], vtype=gp.GRB.INTEGER)
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T:
                    self.x_fqu[q, u, t].setAttr('VarName', f"x_fqu{q}{u}{t}")
        # total amount of fresh HOG 'c' transported from wholesaler 'u' to retailer 'r'
        self.x_cur = self.m.addVars([(u, r, t) for u in self.par.U for r in self.par.R for t in self.par.T] , vtype=gp.GRB.INTEGER)
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.x_cur[u, r, t].setAttr('VarName', f"x_cur{u}{r}{t}")
        # total amount of fillet 'e' transported from wholesaler 'u' to retailer 'r'
        self.x_eur = self.m.addVars([(u, r, t) for u in self.par.U for r in self.par.R for t in self.par.T] , vtype=gp.GRB.INTEGER)
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.x_eur[u, r, t].setAttr('VarName', f"x_eur{u}{r}{t}")
        # total amount of by-product 'f' transported from wholesaler 'u' to retailer 'r'
        self.x_fur = self.m.addVars([(u, r, t) for u in self.par.U for r in self.par.R for t in self.par.T] , vtype=gp.GRB.INTEGER)
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.x_fur[u, r, t].setAttr('VarName', f"x_fur{u}{r}{t}")

    def get_1value(self, df, index, index_val):
        df_raw = df[df[index] == index_val]
        return df_raw[df.columns[-1]].iloc[0] if not df_raw.empty else 0

    def get_2value(self, df, index1, index2, index_val1, index_val2):
        df_raw = df[(df[index1] == index_val1) & (df[index2] == index_val2)]
        return df_raw[df.columns[-1]].iloc[0] if not df_raw.empty else 0
    
    def get_3value(self, df, index1, index2, index3, index_val1, index_val2, index_val3):
        df_raw = df[(df[index1] == index_val1) & (df[index2] == index_val2) & (df[index3] == index_val3)]
        return df_raw[df.columns[-1]].iloc[0] if not df_raw.empty else 0


    def set_obj(self, p):
        # transportation cost  
        self.transportation_cost = gp.quicksum(
            self.get_3value(self.par.G_aijt, 'i', 'j', 't', i, j, t) * self.x_a[i, j, t] for i in self.par.I for j in self.par.J for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_bjpt, 'j', 'p', 't', j, p, t) * self.x_b[j, p, t] for j in self.par.J for p in self.par.P for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_cpqt, 'p', 'q', 't', p, q, t) * self.x_cpq[p, q, t] for p in self.par.P for q in self.par.Q for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_cput, 'p', 'u', 't', p, u, t) * self.x_cpu[p, u, t] for p in self.par.P for u in self.par.U for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_equt, 'q', 'u', 't', q, u, t) * self.x_equ[q, u, t] + self.get_3value(self.par.G_fqut, 'q', 'u', 't', q, u, t) * self.x_fqu[q, u, t] for q in self.par.Q for u in self.par.U for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_curt, 'u', 'r', 't', u, r, t) * self.x_cur[u, r, t] for u in self.par.U for r in self.par.R for t in self.par.T
            ) + gp.quicksum(
                self.get_3value(self.par.G_eurt, 'u', 'r', 't', u, r, t) * self.x_eur[u, r, t] + self.get_3value(self.par.G_furt, 'u', 'r', 't', u, r, t) * self.x_fur[u, r, t] for u in self.par.U for r in self.par.R for t in self.par.T
            )

        # fuel cost
        self.fuel_cost = gp.quicksum(
            self.par.a_t * self.get_2value(self.par.F_aij, 'i', 'j', i, j) * self.get_2value(self.par.W_ij, 'i', 'j', i, j) * self.x_a[i, j, t] for i in self.par.I for j in self.par.J for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * self.get_2value(self.par.F_bjp, 'j', 'p', j, p) * self.get_2value(self.par.W_jp, 'j', 'p', j, p) * self.x_b[j, p, t] for j in self.par.J for p in self.par.P for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * self.get_2value(self.par.F_cpq, 'p', 'q', p, q) * self.get_2value(self.par.W_pq, 'p', 'q', p, q) * self.x_cpq[p, q, t] for p in self.par.P for q in self.par.Q for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * (self.get_2value(self.par.F_equ, 'q', 'u', q, u) * self.x_equ[q, u, t] + self.get_2value(self.par.F_fqu, 'q', 'u', q, u) * self.x_fqu[q, u, t]) * self.get_2value(self.par.W_qu, 'q', 'u', q, u) for q in self.par.Q for u in self.par.U for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * (self.get_2value(self.par.F_eur, 'u', 'r', u, r) * self.x_eur[u, r, t] + self.get_2value(self.par.F_fur, 'u', 'r', u, r) * self.x_fur[u, r, t]) * self.get_2value(self.par.W_ur, 'u', 'r', u, r) for u in self.par.U for r in self.par.R for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * self.get_2value(self.par.F_cpu, 'p', 'u', p, u) * self.get_2value(self.par.W_pu, 'p', 'u', p, u) * self.x_cpu[p, u, t] for p in self.par.P for u in self.par.U for t in self.par.T
            ) + gp.quicksum(
                self.par.a_t * self.get_2value(self.par.F_cur, 'u', 'r', u, r) * self.get_2value(self.par.W_ur, 'u', 'r', u, r) * self.x_cur[u, r, t] for u in self.par.U for r in self.par.R for t in self.par.T
            )

        # inventory cost
        self.inventory_cost = gp.quicksum(
            self.get_2value(self.par.H_bjt, 'j', 't', j, t) * self.ip_b[j,t] for j in self.par.J for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.H_cpt, 'p', 't', p, t) * self.ip_cp[p, t] for p in self.par.P for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.H_cut, 'u','t', u, t) * self.ip_cu[u, t] for u in self.par.U for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.H_eqt, 'q','t', q, t) * self.ip_eq[q, t] + self.get_2value(self.par.H_fqt, 'q','t', q, t) * self.ip_fq[q, t] for q in self.par.Q for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.H_eut, 'u', 't', u, t) * self.ip_eu[u, t] + self.get_2value(self.par.H_fut, 'u','t', u, t) * self.ip_fu[u, t] for u in self.par.U for t in self.par.T
            )

        # processing cost
        self.processing_cost = gp.quicksum(
            self.get_2value(self.par.PC_bjt, 'j', 't', j, t) * self.tp_b[j, t] for j in self.par.J for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.PC_cpt, 'p', 't', p,t) * self.tp_c[p, t] for p in self.par.P for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.PC_eqt, 'q', 't', q, t) * self.tp_e[q, t] + self.get_2value(self.par.PC_fqt, 'q', 't', q, t) * self.tp_f[q, t] for q in self.par.Q for t in self.par.T
            )

        # residual cost
        self.residual_cost = gp.quicksum(
            self.get_2value(self.par.PW_ajt, 'j', 't', j, t) * self.tw_a[j, t] for j in self.par.J for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.PW_bpt, 'p', 't', p, t) * self.tw_b[p, t] for p in self.par.P for t in self.par.T
            ) + gp.quicksum(
                self.get_2value(self.par.PW_cqt, 'q', 't', q, t) * self.tw_c[q, t] for q in self.par.Q for t in self.par.T
            )

        # objective function
        self.m.setObjective(
            self.transportation_cost + self.fuel_cost + self.inventory_cost + self.processing_cost + self.residual_cost,
            gp.GRB.MINIMIZE
        )

    def set_constraints(self, p):
        # (7) the n° of live salmon from each salmon farm shipped to the slaughterhouse must be equal to the capacity of the farm
        for i in self.par.I:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_a[i, j, t] for j in self.par.J) <= self.get_2value(self.par.AC_ait, 'i', 't', i, t), name="Constr7")
        # (8) the total n° of HOG and the residual obtained depends on the total amount of live salmon receveid at the slaughterhouse
        for j in self.par.J:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_a[i, j, t] for i in self.par.I) == (self.tp_b[j, t] + self.tw_a[j, t]), name="Constr8")
        # (9) the total amount of HOG processed after primary proc. pl. depends on the quantity of HOG received from slaughterhouse
        for p in self.par.P:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_b[j, p, t] for j in self.par.J) == (self.tp_c[p, t] + self.tw_b[p, t]), name="Constr9")
        # (10) the total amount of fillets & salmon by-product after secondary proc. pl. depends on the amount of HOG receveid from primary proc. pl.
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_cpq[p, q, t] for p in self.par.P) == (self.tp_e[q, t] + self.tp_f[q, t] + self.tw_c[q, t]), name="Constr10")
        # 3 constraints for residuals:
        # (11) from live to HOG salmon
        for j in self.par.J:
            for t in self.par.T:
                self.m.addConstr(0.05 * gp.quicksum(self.x_a[i, j, t] for i in self.par.I) <= self.tw_a[j, t], name="Constr11a")
                self.m.addConstr(0.2 * gp.quicksum(self.x_a[i, j, t] for i in self.par.I) >= self.tw_a[j, t], name="Constr11b")
        # (12) from HOG to fillets
        for p in self.par.P:
            for t in self.par.T:
                self.m.addConstr(0.05 * gp.quicksum(self.x_b[j, p, t] for j in self.par.J) <= self.tw_b[p, t], name="Constr12a")
                self.m.addConstr(0.2 * gp.quicksum(self.x_b[j, p, t] for j in self.par.J) >= self.tw_b[p, t], name="Constr12b")
        # (13) from HOG to salmon by-products
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(0.05 * gp.quicksum(self.x_cpq[p, q, t] for p in self.par.P) <= self.tw_c[q, t], name="Constr13a")
                self.m.addConstr(0.2 * gp.quicksum(self.x_cpq[p, q, t] for p in self.par.P) >= self.tw_c[q, t], name="Constr13b")
        # (14) the sum of the available inventory of HOG and the total amount of HOG processed should be less or equal to the maximum storage capacity of HOG at slaughterhouse
        for j in self.par.J:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.tp_b[j, t] <= self.get_2value(self.par.Cap_bjt, 'j', 't', j, t) - self.ip_b[j,t-1])
                else:
                    #if self.ip_b[j,t] == 0:
                    self.m.addConstr(self.tp_b[j, t] <= self.get_2value(self.par.Cap_bjt, 'j', 't', j, t))
          

        # (15) the sum of the available inventory of HOG and the total amount of HOG processed shold be equal or less then the maximum storage capacity for HOG at the primary proc. pl.
        for p in self.par.P:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.tp_c[p, t] <= self.get_2value(self.par.Cap_cpt, 'p', 't', p, t) - self.ip_cp[p,t-1], name="Constr15b")
                else:
                    #if self.ip_cp[p,t] == 0:
                    self.m.addConstr(self.tp_c[p, t] <= self.get_2value(self.par.Cap_cpt, 'p', 't', p, t), name="Constr15a")
                    
        # (16-17) " " for both fillets & salmon by-products
        for q in self.par.Q:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.tp_e[q, t] <= self.get_2value(self.par.Cap_eqt, 'q', 't', q, t)  - self.ip_eq[q,t-1], name="Constr16b")
                    self.m.addConstr(self.tp_f[q, t] <= self.get_2value(self.par.Cap_fqt, 'q', 't', q, t) - self.ip_fq[q,t-1], name="Constr17b")
                else:
                    #if self.ip_eq[q,t] == 0:
                    self.m.addConstr(self.tp_e[q, t] <= self.get_2value(self.par.Cap_eqt, 'q', 't', q, t), name="Constr16a")
                    #elif self.ip_fq[q,t] == 0:    
                    self.m.addConstr(self.tp_f[q, t] <= self.get_2value(self.par.Cap_fqt, 'q', 't', q, t), name="Constr17a")
             
                    
        # (18) the sum of the total HOG from primary to a wholesaler and the available inventory of HOG at the wholesaler whould be less or equal to the maximum storage capacity of HOG at the wholesaler
        for u in self.par.U:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(gp.quicksum(self.x_cpu[p, u, t] for p in self.par.P) <= self.get_2value(self.par.Cap_cut, 'u', 't', u, t)  - self.ip_cu[u,t-1], name="Constr18b")
                else:
                    #if self.ip_cu[u,t] == 0:
                    self.m.addConstr(gp.quicksum(self.x_cpu[p, u, t] for p in self.par.P) <= self.get_2value(self.par.Cap_cut, 'u', 't', u, t), name="Constr18a")
                    
        # (19-20) storage capacity of the wholesaler for fillets & salmon by-products
        for u in self.par.U:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(gp.quicksum(self.x_equ[q, u, t] for q in self.par.Q) <= self.get_2value(self.par.Cap_eut, 'u', 't', u, t) - self.ip_eu[u,t-1], name="Constr19")
                    self.m.addConstr(gp.quicksum(self.x_fqu[q, u, t] for q in self.par.Q) <= self.get_2value(self.par.Cap_fut, 'u', 't', u, t) - self.ip_fu[u,t-1], name="Constr20")
                else:
                    #if self.ip_eu[u,t] == 0:
                    self.m.addConstr(gp.quicksum(self.x_equ[q, u, t] for q in self.par.Q) <= self.get_2value(self.par.Cap_eut, 'u', 't', u, t), name="Constr19a")
                    #el#ip_fu[u,t] == 0:
                    self.m.addConstr(gp.quicksum(self.x_fqu[q, u, t] for q in self.par.Q) <= self.get_2value(self.par.Cap_fut, 'u', 't', u, t), name="Constr20a")
         
                    
        # (21) inventory balancing constraint for HOG salmon at the slaughterhouse
        for j in self.par.J:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.ip_b[j, t] == (self.tp_b[j, t] - gp.quicksum(self.x_b[j, p, t] for p in self.par.P)) + self.ip_b[j,t-1], name="Constr21b")
                else:
                    #ip_b[j,t] == 0:
                    self.m.addConstr(self.ip_b[j, t] == (self.tp_b[j, t] - gp.quicksum(self.x_b[j, p, t] for p in self.par.P)), name="Constr21a") 
        
        # (22-23) inventory balancing constraints for fresh HOG salmon products at the primary processing plant and the wholesaler
        for p in self.par.P:
            for t in self.par.T:
                if t > 1: 
                    self.m.addConstr(self.ip_cp[p, t] == (self.tp_c[p, t] - gp.quicksum(self.x_cpq[p, q, t] for q in self.par.Q) - gp.quicksum(self.x_cpu[p, u, t] for u in self.par.U)) + self.ip_cp[p,t-1], name="Constr22a")
                else:
                    #if self.ip_cp[p,t] == 0:
                    self.m.addConstr(self.ip_cp[p, t] == (self.tp_c[p, t] - gp.quicksum(self.x_cpq[p, q, t] for q in self.par.Q) - gp.quicksum(self.x_cpu[p, u, t] for u in self.par.U)), name="Constr22b") 
                
        for t in self.par.T:
            for u in self.par.U:
                if t > 1:
                    self.m.addConstr(self.ip_cu[u, t] == (gp.quicksum(self.x_cpu[p, u, t] for p in self.par.P) - gp.quicksum(self.x_cur[u, r, t] for r in self.par.R) + self.ip_cu[u,t-1]), name="Constr23b")
                else: 
                     #if self.ip_cu[u,t] == 0:
                    self.m.addConstr(self.ip_cu[u, t] == (gp.quicksum(self.x_cpu[p, u, t] for p in self.par.P) - gp.quicksum(self.x_cur[u, r, t] for r in self.par.R)), name="Constr23a")
                    
        # (24-25) inventory balancing constraints at the secondary processing plant for whole fillets and salmon by-products
        for q in self.par.Q:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.ip_eq[q, t] == (self.tp_e[q, t] - gp.quicksum(self.x_equ[q, u, t] for u in self.par.U) + self.ip_eq[q,t-1]), name="Constr24a")
                    self.m.addConstr(self.ip_fq[q, t] == (self.tp_f[q, t] - gp.quicksum(self.x_fqu[q, u, t] for u in self.par.U)+ self.ip_fq[q,t-1]), name="Constr25a")
                else:
                    #if self.ip_eq[q,t] == 0:
                    self.m.addConstr(self.ip_eq[q, t] == (self.tp_e[q, t] - gp.quicksum(self.x_equ[q, u, t] for u in self.par.U)), name="Constr24b")
                    self.m.addConstr(self.ip_fq[q, t] == (self.tp_f[q, t] - gp.quicksum(self.x_fqu[q, u, t] for u in self.par.U)), name="Constr25b")
                    
         # (26-27) inventory balancing constraints for whole fillets and salmon by-products respectively at the wholesaler
        for u in self.par.U:
            for t in self.par.T:
                if t > 1:
                    self.m.addConstr(self.ip_eu[u, t] == (gp.quicksum(self.x_equ[q, u, t] for q in self.par.Q) - gp.quicksum(self.x_eur[u, r, t] for r in self.par.R) + self.ip_eu[u, t-1]), name="Constr26b")
                    self.m.addConstr(self.ip_fu[u, t] == (gp.quicksum(self.x_fqu[q, u, t] for q in self.par.Q) - gp.quicksum(self.x_fur[u, r, t] for r in self.par.R) + self.ip_fu[u, t-1]), name="Constr27b")

                else:
                    self.m.addConstr(self.ip_eu[u, t] == (gp.quicksum(self.x_equ[q, u, t] for q in self.par.Q) - gp.quicksum(self.x_eur[u, r, t] for r in self.par.R)), name="Constr26a")
                    self.m.addConstr(self.ip_fu[u, t] == (gp.quicksum(self.x_fqu[q, u, t] for q in self.par.Q) - gp.quicksum(self.x_fur[u, r, t] for r in self.par.R)), name="Constr27a")
                
                    
        
        # (28) demand constraints at the retail level for fresh HOG
        for r in self.par.R:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_cur[u, r, t] for u in self.par.U) == self.get_2value(self.par.D_crt, 'r', 't', r, t), name="Constr28")
        # (29) demand constraints at the retail level for fillets
        for r in self.par.R:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_eur[u, r, t] for u in self.par.U) == self.get_2value(self.par.D_ert, 'r', 't', r, t), name="Constr29")
        # (30) demand constraints at the retail level for salmon by-products
        for r in self.par.R:
            for t in self.par.T:
                self.m.addConstr(gp.quicksum(self.x_fur[u, r, t] for u in self.par.U) == self.get_2value(self.par.D_frt, 'r', 't', r, t), name="Constr30")

        # (31) the overall carbon emissions emitted from transportation of various salmon products should be less than or equal to the maximum allowable carbon emission limit
        self.carbon_emission = self.par.E_CO2 * (gp.quicksum(
                self.get_2value(self.par.F_aij, 'i', 'j', i, j) * self.get_2value(self.par.W_ij, 'i', 'j', i, j) * self.x_a[i, j, t] for i in self.par.I for j in self.par.J
                ) + gp.quicksum(
                    self.get_2value(self.par.F_bjp, 'j', 'p', j, p) * self.get_2value(self.par.W_jp, 'j', 'p', j, p) * self.x_b[j, p, t] for j in self.par.J for p in self.par.P
                ) + gp.quicksum(
                    self.get_2value(self.par.F_cpq, 'p', 'q', p, q) * self.get_2value(self.par.W_pq, 'p', 'q', p, q) * self.x_cpq[p, q, t] for p in self.par.P for q in self.par.Q
                ) + gp.quicksum(
                    (self.get_2value(self.par.F_equ, 'q', 'u', q, u) * self.x_equ[q, u, t] + self.get_2value(self.par.F_fqu, 'q', 'u', q, u) * self.x_fqu[q, u, t]) * self.get_2value(self.par.W_qu, 'q', 'u', q, u) for q in self.par.Q for u in self.par.U
                ) + gp.quicksum(
                    (self.get_2value(self.par.F_eur, 'u', 'r', u, r) * self.x_eur[u, r, t] + self.get_2value(self.par.F_fur, 'u', 'r', u, r) * self.x_fur[u, r, t]) * self.get_2value(self.par.W_ur, 'u', 'r', u, r) for u in self.par.U for r in self.par.R
                ) + gp.quicksum(
                    self.get_2value(self.par.F_cpu, 'p', 'u', p, u) * self.get_2value(self.par.W_pu, 'p', 'u', p, u) * self.x_cpu[p, u, t] for p in self.par.P for u in self.par.U
                ) + gp.quicksum(
                    self.get_2value(self.par.F_cur, 'u', 'r', u, r) * self.get_2value(self.par.W_ur, 'u', 'r', u, r) * self.x_cur[u, r, t] for u in self.par.U for r in self.par.R
                )
            )
        self.m.addConstr(self.carbon_emission <= self.par.E_max, name="Constr31"
        )


        # (B1) the number of products flowing at each transportation link from farms to slaughterhouses should be less than or equal to the maximum transportation capacity at each transportation link
        for i in self.par.I:
            for j in self.par.J:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_aijt, 'i', 'j', 't', i, j, t) == 0:
                        self.m.addConstr(self.x_a[i, j, t] == 0)
                    else:
                        self.m.addConstr(self.x_a[i, j, t] <= self.get_3value(self.par.TC_aijt, 'i', 'j', 't', i, j, t))
        
        # (B2) the number of products flowing at each transportation link from slaughterhouses to primary processing plants should be less than or equal to the maximum transportation capacity at each transportation link
        for j in self.par.J:
            for p in self.par.P:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_bjpt, 'j', 'p', 't', j, p, t) == 0:
                        self.m.addConstr(self.x_b[j, p, t] == 0)
                    else:
                        self.m.addConstr(self.x_b[j, p, t] <= self.get_3value(self.par.TC_bjpt, 'j', 'p', 't', j, p, t))     

        # (B3) the number of fresh HOG salmon flowing at each transportation link from primary to secondary processing plants should be less than or equal to the maximum transportation capacity at each transportation link
        for p in self.par.P:
            for q in self.par.Q:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_cpqt, 'p', 'q', 't', p, q, t) == 0:
                        self.m.addConstr(self.x_cpq[p, q, t] == 0)
                    else:
                        self.m.addConstr(self.x_cpq[p, q, t] <= self.get_3value(self.par.TC_cpqt, 'p', 'q', 't', p, q, t))    

        # (B4) the number of fresh HOG salmon flowing at each transportation link from primary processing plants to wholesalers should be less than or equal to the maximum transportation capacity at each transportation link
        for p in self.par.P:
            for u in self.par.U:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_cput, 'p', 'u', 't', p, u, t) == 0:
                        self.m.addConstr(self.x_cpu[p, u, t] == 0)
                    else:
                        self.m.addConstr(self.x_cpu[p, u, t] <= self.get_3value(self.par.TC_cput, 'p', 'u', 't', p, u, t))

        # (B5) the number of whole fillets flowing at each transportation link from secondary processing plants to wholesalers should be less than or equal to the maximum transportation capacity at each transportation link
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_equt, 'q', 'u', 't', q, u, t) == 0:
                        self.m.addConstr(self.x_equ[q, u, t] == 0)
                    else:
                        self.m.addConstr(self.x_equ[q, u, t] <= self.get_3value(self.par.TC_equt, 'q', 'u', 't', q, u, t))  


        # (B6) the number of salmon by-products flowing at each transportation link from secondary processing plants to wholesalers should be less than or equal to the maximum transportation capacity at each transportation link
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_fqut, 'q', 'u', 't', q, u, t) == 0:
                        self.m.addConstr(self.x_fqu[q, u, t] == 0)
                    else:
                        self.m.addConstr(self.x_fqu[q, u, t] <= self.get_3value(self.par.TC_fqut, 'q', 'u', 't', q, u, t))  

        # (B7) the number of fresh HOG products flowing at each transportation link from wholesalers to retailers should be less than or equal to the maximum transportation capacity at each transportation link
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_curt, 'u', 'r', 't', u, r, t) == 0:
                        self.m.addConstr(self.x_cur[u, r, t] == 0, name="ConstrB7")
                    else:
                        self.m.addConstr(self.x_cur[u, r, t] <= self.get_3value(self.par.TC_curt, 'u', 'r', 't', u, r, t))  


        # (B8) the number of fresh whole fillets flowing at each transportation link from wholesalers to retailers should be less than or equal to the maximum transportation capacity at each transportation link
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_eurt, 'u', 'r', 't', u, r, t) == 0:
                        self.m.addConstr(self.x_eur[u, r, t] == 0)
                    else:
                        self.m.addConstr(self.x_eur[u, r, t] <= self.get_3value(self.par.TC_eurt, 'u', 'r', 't', u, r, t))  


        # (B9) the number of fresh salmon by-products flowing at each transportation link from wholesalers to retailers should be less than or equal to the maximum transportation capacity at each transportation link
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T: 
                    if  self.get_3value(self.par.TC_furt, 'u', 'r', 't', u, r, t) == 0:
                        self.m.addConstr(self.x_fur[u, r, t] == 0, name="ConstrB9")
                    else:
                        self.m.addConstr(self.x_fur[u, r, t] <= self.get_3value(self.par.TC_furt, 'u', 'r', 't', u, r, t))  


        # (B10a) total processed amount of HOG product 'b' at slaughterhouse j in period t is equal or larger than 0
        for j in self.par.J:
            for t in self.par.T:
                self.m.addConstr(self.tp_b[j, t] >= 0)


        # (B10b) total processed amount of HOG product 'c' at primary processing plant p in period t is equal or larger than 0
        for p in self.par.P:
            for t in self.par.T:
                self.m.addConstr(self.tp_c[p, t] >= 0)

        # (B10c) total processed amount of whole fillets 'e' at secondary processing plant q in period t is equal or larger than 0
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(self.tp_e[q, t] >= 0)

        # (B10d) total processed amount salmon by-products 'f' at secondary processing plant q in period t is equal or larger than 0
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(self.tp_f[q, t] >= 0)

        # (B10e) total amount of residual obtained after processing live salmon product 'a' at slaughterhouse j in period t is equal or larger than 0
        for j in self.par.J:
            for t in self.par.T:
                self.m.addConstr(self.tw_a[j, t] >= 0)

        # (B10f) total amount of residual obtained after processing HOG product 'b' at primary processing p in period t is equal or larger than 0 
        for p in self.par.P:
            for t in self.par.T:
                self.m.addConstr(self.tw_b[p, t] >= 0)
        
        # (B10g) total amount of residual obtained after processing HOG product 'c' at secondary processing plant q in period t is equal or larger than 0
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(self.tw_c[q, t] >= 0)


        # (B11a) inventory of HOG product 'b' at slaughterhouse j in period t is equal or larger than 0
        for j in self.par.J:
            for t in self.par.T:
                self.m.addConstr(self.ip_b[j, t] >= 0)

        # (B11b) inventory of HOG product 'c' at primary processing plant p in period t is equal or larger than 0
        for p in self.par.P:
            for t in self.par.T:
                self.m.addConstr(self.ip_cp[p, t] >= 0)

        # (B11c) inventory of HOG product 'c' at wholesaler u in period t is equal or larger than 0
        for u in self.par.U:
            for t in self.par.T:
                self.m.addConstr(self.ip_cu[u, t] >= 0)

        # (B11d) inventory of whole fillets 'e' at secondary processing plant q in period t is equal or larger than 0
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(self.ip_eq[q, t] >= 0)

        # (B11e) inventory of salmon by-products 'f' at secondary processing plant q in period t is equal or larger than 0
        for q in self.par.Q:
            for t in self.par.T:
                self.m.addConstr(self.ip_fq[q, t] >= 0)

        # (B11f) inventory of whole fillets 'e' at wholesaler u in period t is equal or larger than 0 
        for u in self.par.U:
            for t in self.par.T:
                self.m.addConstr(self.ip_eu[u, t] >= 0)
        
        # (B11g) inventory of salmon by-producst 'f' at wholesaler u in period t is equal or larger than 0 
        for u in self.par.U:
            for t in self.par.T:
                self.m.addConstr(self.ip_fu[u, t] >= 0)



        # (B12a) total amount of live salmon 'a' transported from salmon farm i to slaughterhouse j is equal or larger than 0
        for i in self.par.I:
            for j in self.par.J:
                for t in self.par.T:
                    self.m.addConstr(self.x_a[i, j, t] >= 0)


        # (B12b) total amount of HOG product 'b' transported from slaughterhouse j to primary processing plant p is equal or larger than 0
        for j in self.par.J:
            for p in self.par.P:
                for t in self.par.T:
                    self.m.addConstr(self.x_b[j, p, t] >= 0)

        # (B12c) total amount of HOG product 'c' transported from primary processing plant p to secondary processing plant q is equal or larger than 0
        for p in self.par.P:
            for q in self.par.Q:
                for t in self.par.T:
                    self.m.addConstr(self.x_cpq[p, q, t] >= 0)


        # (B12d) total amount of HOG product 'c' transported from primary processing plant p to wholesaler u is equal or larger than 0
        for p in self.par.P:
            for u in self.par.U:
                for t in self.par.T:
                    self.m.addConstr(self.x_cpu[p, u, t] >= 0)

        # (B12e) total amount of whole fillets 'e' transported from secondary processing plant q to wholesaler u is equal or larger than 0
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T:
                    self.m.addConstr(self.x_equ[q, u, t] >= 0)

        # (B12f) total amount of salmon by-products 'f' transported from secondary processing plant q to wholesaler u is equal or larger than 0
        for q in self.par.Q:
            for u in self.par.U:
                for t in self.par.T:
                    self.m.addConstr(self.x_fqu[q, u, t] >= 0)

        # (B12g) total amount of HOG product 'c' transported from wholesaler u to retailer r is equal or larger than 0
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.m.addConstr(self.x_cur[u, r, t] >= 0)

        # (B12h) total amount of whole fillets 'e' transported from wholesaler u to retailer r is equal or larger than 0
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.m.addConstr(self.x_eur[u, r, t] >= 0)
        
        # (B12i) total amount of salmon by-products 'f' transported from wholesaler u to retailer r is equal or larger than 0
        for u in self.par.U:
            for r in self.par.R:
                for t in self.par.T:
                    self.m.addConstr(self.x_fur[u, r, t] >= 0)
  

    def optimise(self, logs=False):
        if not logs:
            self.m.setParam('OutputFlag', 0)
        self.m.optimize()
        print(self.m.status)

        if self.m.status != gp.GRB.INFEASIBLE:
            # for var in self.m.getVars():
            #     print(var.varName, var.x)
            print('\nNumber of decision variables:', self.m.NumVars)
            print('Number of constraints:', self.m.NumConstrs)
            print('Size of the input data:', self.par.num_data_points, 'data points')
            print('Objective function:', self.m.objVal)
            print('\nFuel cost:', self.fuel_cost.getValue())
            print('Residual cost:', self.residual_cost.getValue())
            print('Inventory cost:', self.inventory_cost.getValue())
            print('Transportation cost:', self.transportation_cost.getValue())
            print('Processing cost:', self.processing_cost.getValue())
            print('Carbon emission:', self.carbon_emission.getValue(), 'KgCO2')


        if self.m.status == gp.GRB.INFEASIBLE:
            self.m.computeIIS()
            print('The following constraints and variables are in the IIS:')
            for c in self.m.getConstrs():
                if c.IISConstr: print(f'\t{c.constrname}: {self.m.getRow(c)} {c.Sense} {c.RHS}')

        