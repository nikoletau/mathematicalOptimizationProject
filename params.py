import pandas as pd

class Params:
    def __init__(self):

        # PATH FROM FARM 'i' TO SLAUGHTERHOUSE 'j'
        ij = pd.read_csv("istanza1/ij.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of Live Salmon 'a'
        self.TC_aijt = ij.loc[:, ['i', 'j', 't', 'TC_aijt']].drop_duplicates()
        # transportation cost for the shipment of per unit Live Salmon
        self.G_aijt = ij.loc[:, ['i', 'j', 't', 'G_aijt']].drop_duplicates()
        # fuel consumed in shipping Live Salmon 'a' via certain mode of transport
        self.F_aij = ij.loc[:, ['i', 'j', 't', 'F_aij']].drop_duplicates()
        # distance
        self.W_ij = ij.loc[:, ['i', 'j', 'W_ij']].drop_duplicates()

        # PATH FROM SLAUGHTERHOUSE 'j' TO PRIMARY PROCESSING PLANT 'p'
        jp = pd.read_csv("istanza1/jp.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of HOG 'b'
        self.TC_bjpt = jp.loc[:, ['j', 'p', 't', 'TC_bjpt']].drop_duplicates()
        # transportation cost for the shipment of per unit HOG 'b'
        self.G_bjpt = jp.loc[:, ['j', 'p', 't', 'G_bjpt']].drop_duplicates()
        # fuel consumed in shipping HOG 'b'
        self.F_bjp = jp.loc[:, ['j', 'p', 'F_bjp']].drop_duplicates()
        # distance
        self.W_jp = jp.loc[:, ['j', 'p', 'W_jp']].drop_duplicates()

        # PATH FROM PRIMARY PROCESSING PLANT 'p' TO SECONDARY PROCESSING PLANT 'q'
        pq = pd.read_csv("istanza1/pq.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of fresh HOG 'c'
        self.TC_cpqt = pq.loc[:, ['p', 'q', 't', 'TC_cpqt']].drop_duplicates()
        # transportation cost for the shipment of per unit fresh HOG 'c'
        self.G_cpqt = pq.loc[:, ['p', 'q', 't', 'G_cpqt']].drop_duplicates()
        # fuel consumed in shipping fresh HOG 'c'
        self.F_cpq = pq.loc[:, ['p', 'q', 'F_cpq']].drop_duplicates()
        # distance
        self.W_pq = pq.loc[:, ['p', 'q', 'W_pq']].drop_duplicates()

        # PATH FROM PRIMARY PROCESSING PLANT 'p' TO WHOLESALER 'u'
        pu = pd.read_csv("istanza1/pu.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of fresh HOG 'c'
        self.TC_cput = pu.loc[:, ['p', 'u', 't', 'TC_cput']].drop_duplicates()
        # transportation cost for the shipment of per unit fresh HOG 'c'
        self.G_cput = pu.loc[:, ['p', 'u', 't', 'G_cput']].drop_duplicates()
        # fuel consumed in shipping fresh HOG 'c'
        self.F_cpu = pu.loc[:, ['p', 'u', 'F_cpu']].drop_duplicates()
        # distance
        self.W_pu = pu.loc[:, ['p', 'u', 'W_pu']].drop_duplicates()

        # PATH FROM SECONDARY PROCESSING PLANT 'q' TO wholesaler 'u'
        qu = pd.read_csv("istanza1/qu.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of fillet 'e'
        self.TC_equt = qu.loc[:, ['q', 'u', 't', 'TC_equt']].drop_duplicates()
        # capacity of the transportation mode for the shipment of by-product 'f'
        self.TC_fqut = qu.loc[:, ['q', 'u', 't', 'TC_fqut']].drop_duplicates()
        # transportation cost for the shipment of per unit fillet 'e'
        self.G_equt = qu.loc[:, ['q', 'u', 't', 'G_equt']].drop_duplicates()
        # transportation cost for the shipment of per unit by-product 'f'
        self.G_fqut = qu.loc[:, ['q', 'u', 't', 'G_fqut']].drop_duplicates()
        # fuel consumed in shipping fillet 'e'
        self.F_equ = qu.loc[:, ['q', 'u', 'F_equ']].drop_duplicates()
        # fuel consumed in shipping by-product 'f'
        self.F_fqu = qu.loc[:, ['q', 'u', 'F_fqu']].drop_duplicates()
        # distance
        self.W_qu = qu.loc[:, ['q', 'u', 'W_qu']].drop_duplicates()

        # PATH FROM WHOLESALER 'u' TO RETAILER 'r'
        ur = pd.read_csv("istanza1/ur.csv", sep=";", header=0)
        # capacity of the transportation mode for the shipment of fresh HOG 'c'
        self.TC_curt = ur.loc[:, ['u', 'r', 't', 'TC_curt']].drop_duplicates()
        # capacity of the transportation mode for the shipment of fillet 'e'
        self.TC_eurt = ur.loc[:, ['u', 'r', 't', 'TC_eurt']].drop_duplicates()
        # capacity of the transportation mode for the shipment of by-product 'f'
        self.TC_furt = ur.loc[:, ['u', 'r', 't', 'TC_furt']].drop_duplicates()
        # transportation cost for the shipment of per unit fresh HOG 'c'
        self.G_curt = ur.loc[:, ['u', 'r', 't', 'G_curt']].drop_duplicates()
        # transportation cost for the shipment of per unit fillet 'e'
        self.G_eurt = ur.loc[:, ['u', 'r', 't', 'G_eurt']].drop_duplicates()
        # transportation cost for the shipment of per unit by-product 'f'
        self.G_furt = ur.loc[:, ['u', 'r', 't', 'G_furt']].drop_duplicates()
        # fuel consumed in shipping fresh HOG 'c'
        self.F_cur = ur.loc[:, ['u', 'r', 'F_cur']].drop_duplicates()
        # fuel consumed in shipping fillet 'e'
        self.F_eur = ur.loc[:, ['u', 'r', 'F_eur']].drop_duplicates()
        # fuel consumed in shipping by-product 'f'
        self.F_fur = ur.loc[:, ['u', 'r', 'F_fur']].drop_duplicates()
        # distance
        self.W_ur = ur.loc[:, ['u', 'r', 'W_ur']].drop_duplicates()

        # FARM 'i'
        farm = pd.read_csv("istanza1/farm.csv", sep=";", header=0)
        # capacity of Live Salmon 'a'
        self.AC_ait = farm.loc[:, ['i', 't', 'AC_ait']].drop_duplicates()

        # SLAUGHTERHOUSE 'j'
        sl = pd.read_csv("istanza1/slaughterhouse.csv", sep=";", header=0)
        # maximum storage capacity of HOG 'b'
        self.Cap_bjt = sl.loc[:, ['j', 't', 'Cap_bjt']].drop_duplicates()
        # inventory holding cost per unit of HOG 'b'
        self.H_bjt = sl.loc[:, ['j', 't', 'H_bjt']].drop_duplicates()
        # processing cost for per unit of HOG 'b'
        self.PC_bjt = sl.loc[:, ['j', 't', 'PC_bjt']].drop_duplicates()
        # residual cost for per unit of residual amount obtained after processing Live Salmon 'a'
        self.PW_ajt = sl.loc[:, ['j', 't', 'PW_ajt']].drop_duplicates()

        # PRIMARY PROCESSING PLANT 'p'
        prim = pd.read_csv("istanza1/primary.csv", sep=";", header=0)
        # maximum storage capacity of fresh HOG 'c'
        self.Cap_cpt = prim.loc[:, ['p', 't', 'Cap_cpt']].drop_duplicates()
        # inventory holding cost per unit of fresh HOG 'c'
        self.H_cpt = prim.loc[:, ['p', 't', 'H_cpt']].drop_duplicates()
        # processing cost for per unit of fresh HOG 'c'
        self.PC_cpt = prim.loc[:, ['p', 't', 'PC_cpt']].drop_duplicates()
        # residual cost for per unit of residual amount obtained after processing HOG 'b'
        self.PW_bpt = prim.loc[:, ['p', 't', 'PW_bpt']].drop_duplicates()

        # SECONDARY PROCESSING PLANT 'q'
        sec = pd.read_csv("istanza1/secondary.csv", sep=";", header=0)
        # maximum storage capacity of fillet 'e'
        self.Cap_eqt = sec.loc[:, ['q', 't', 'Cap_eqt']].drop_duplicates()
        # maximum storage capacity of by-product 'f'
        self.Cap_fqt = sec.loc[:, ['q', 't', 'Cap_fqt']].drop_duplicates()
        # inventory holding cost per unit of fillet 'e'
        self.H_eqt = sec.loc[:, ['q', 't', 'H_eqt']].drop_duplicates()
        # maximum storage capacity of by-product 'f'
        self.H_fqt = sec.loc[:, ['q', 't', 'H_fqt']].drop_duplicates()
        # processing cost for per unit of fillet 'e'
        self.PC_eqt = sec.loc[:, ['q', 't', 'PC_eqt']].drop_duplicates()
        # processing cost for per unit of by-product 'f'
        self.PC_fqt = sec.loc[:, ['q', 't', 'PC_fqt']].drop_duplicates()
        # residual cost for per unit of residual amount obtained after processing fresh HOG 'c'
        self.PW_cqt = sec.loc[:, ['q', 't', 'PW_cqt']].drop_duplicates()

        # WHOLESALER 'u'
        whol = pd.read_csv("istanza1/wholesaler.csv", sep=";", header=0)
        # maximum storage capacity of fresh HOG 'c'
        self.Cap_cut = whol.loc[:, ['u', 't', 'Cap_cut']].drop_duplicates()
        # maximum storage capacity of fillet 'e'
        self.Cap_eut = whol.loc[:, ['u', 't', 'Cap_eut']].drop_duplicates()
        # maximum storage capacity of by-product 'f'
        self.Cap_fut = whol.loc[:, ['u', 't', 'Cap_fut']].drop_duplicates()
        # inventory holding cost per unit of fresh HOG 'c'
        self.H_cut = whol.loc[:, ['u', 't', 'H_cut']].drop_duplicates()
        # inventory holding cost per unit of fillet 'e'
        self.H_eut = whol.loc[:, ['u', 't', 'H_eut']].drop_duplicates()
        # inventory holding cost per unit of by-product 'f'
        self.H_fut = whol.loc[:, ['u', 't', 'H_fut']].drop_duplicates()

        # RETAILER 'r'
        ret = pd.read_csv("istanza1/retailer.csv", sep=";", header=0)
        # demand for fresh HOG 'c'
        self.D_crt = ret.loc[:, ['r', 't', 'D_crt']].drop_duplicates()
        # demand for fillet 'e'
        self.D_ert = ret.loc[:, ['r', 't', 'D_ert']].drop_duplicates()
        # demand for by-product 'f'
        self.D_frt = ret.loc[:, ['r', 't', 'D_frt']].drop_duplicates()

        # FIXED PARAMETERS
        self.a_t = 1.5 # fuel price (Euro per litre)
        self.E_CO2 = 2.392 # carbon emission coefficient associated with the fuel
        self.E_max = 1000000 # maximum allowable carbon emission limit

        # RETRIEVE INDICES
        self.I=farm['i'].unique().tolist()
        self.J=sl['j'].unique().tolist()
        self.P=prim['p'].unique().tolist()
        self.Q=sec['q'].unique().tolist()
        self.U=whol['u'].unique().tolist()
        self.R=ret['r'].unique().tolist()
        self.T=farm['t'].unique().tolist()

        # Combine all dataframes in one dataframe to retrieve size of input data
        combined_df = pd.concat([self.AC_ait, self.Cap_bjt, self.Cap_cpt, self.Cap_cut, self.Cap_eqt, self.Cap_eqt, self.Cap_eut, self.Cap_fqt, self.Cap_fut, 
        self.D_crt, self.D_ert, self.D_frt, 
        self.F_aij, self.F_bjp, self.F_cpq, self.F_cpu, self.F_cur, self.F_equ, self.F_eur, self.F_fqu, self.F_fur,
        self.G_aijt, self.G_bjpt, self.G_cpqt, self.G_cput, self.G_curt, self.G_equt, self.G_eurt, self.G_fqut, self.G_furt,
        self.H_bjt, self.H_cpt, self.H_cut, self.H_eqt, self.H_eut, self.H_fqt, self.H_fut,
        self.PC_bjt, self.PC_cpt, self.PC_eqt, self.PC_fqt,
        self.PW_ajt, self.PW_bpt, self.PW_cqt,
        self.TC_aijt, self.TC_bjpt, self.TC_cpqt, self.TC_cput, self.TC_curt, self.TC_equt, self.TC_eurt, self.TC_fqut, self.TC_furt,
        self.W_ij, self.W_jp, self.W_pq, self.W_pu, self.W_qu, self.W_ur])
        self.num_data_points = combined_df.shape[0]
