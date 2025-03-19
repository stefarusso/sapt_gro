import sys
import os
import pandas
#python shift.py final.gro steps
# in the directory itp top mdp files are required

class coord():
    def __init__(self):
        if len(sys.argv) > 2:
            self.n = int(sys.argv[2])
            #self.step = int(sys.argv[3])
            #self.s = self.s/self.step
            self.filename_gro = sys.argv[1]
            self.atom1, self.atom2 = [],[]
            self.x, self.y, self.z = [],[],[]
            self.res_number = []
            self.blank = []
            self.atoms = []
            self.res_name = []
            self.atom_name = []
            self.atom_n = []
            self.cell =[]
            self.deltax, self.deltay, self.deltaz, self.norm = [],[],[],[]
            self.data = []
            self.load()
            self.save_all()
            self.data = pandas.DataFrame(self.data, columns=["Distance","Total_opls","Coulomb_opls","LJ_opls","Total_sapt","Coulomb_sapt","LJ_sapt"] )
            self.analysis()
        else:
            print("Format Required : python script.py file.gro N_steps/2")
    def load(self):
        with open(self.filename_gro) as file:
            blank = file.readline()
            atoms = file.readline()
            atoms = int(atoms.strip().split()[0])
            res_name, atom_name,atom_n,x,y,z = [],[],[],[],[],[]
            res_number = []
            last_resname = []
            count = 0
            res_names = []
            for i in range(atoms):
                line = file.readline()
                res_name_number = line.strip().split()[0]
                res_name.append(line[0:10])
                line=line[11:]
                line=line.strip().split()
                #check if new res_name is found
                if last_resname != res_name_number:
                    #new residue
                    if last_resname:
                        #save count after check is not empty
                        res_number.append(count)
                        res_names.append(last_resname)
                    last_resname = res_name_number
                    count = 1
                elif i==atoms-1:
                    count +=1
                    res_number.append(count)
                    res_names.append(last_resname)
                else:
                    count += 1
                atom_name.append(line[0])
                atom_n.append( int(line[1]) )
                x.append( float(line[2]) )
                y.append( float(line[3]) )
                z.append( float(line[4]) )
            cell = file.readline()
            string1='Residue1:'
            string2='Residue2:'
            print(f'{string1:<12s}{string2:<12s}')
            string1 = "[1-" + str( res_number[0] ) +"]"
            string2 = "["+ str(res_number[0]+1) +"-"+ str(res_number[0]+res_number[1])+"]"
            print(f"{string1:<12}{string2:<12s}")
            while True:
                try:
                    atom1 = int(input("Select Atom Number for Residue1: " ))
                    if atom1 < 1 or atom1 > res_number[0] :
                        print(f"ERROR!!!!  Atom Chosen outside of Residue Range:{atom1} [0-{res_number[0]}]")
                        raise ValueError()
                    else:
                        atom1 += -1
                        break
                except ValueError:
                    pass
            while True:
                try:
                    atom2 = int(input("Select Atom Number for Residue2: " ))
                    if atom2 < res_number[0]+1 or atom2 > res_number[0]+res_number[1] :
                        print(f'ERROR!!!!  Atom Chosen outside of Residue Range:{atom2} [{res_number[0]+1}-{res_number[1]+res_number[0]}]')
                        raise ValueError()
                    else:
                        atom2 += -1
                        break
                except ValueError:
                    pass
            deltax=x[atom2]-x[atom1]
            deltay=y[atom2]-y[atom1]
            deltaz=z[atom2]-z[atom1]
            norm=(deltax**2+deltay**2+deltaz**2)**0.5
            deltax=deltax/norm
            deltay=deltay/norm
            deltaz=deltaz/norm
            self.deltax, self.deltay, self.deltaz, self.norm = deltax, deltay, deltaz, norm
            self.atom1 = atom1
            self.atom2 = atom2
            self.x, self.y, self.z = x, y, z
            self.res_number = res_number
            self.res_names = res_names
            self.res_name = res_name
            self.atom_n = atom_n
            self.atom_name = atom_name
            self.blank = blank
            self.atoms =atoms
            self.cell = cell
    def ask_charge(self):
        q1 = input("Charge of The Residue1 :[1] ").strip()
        m1 = input("Multiplicity of The Residue1 :[1] ").strip()
        q2 = input("Charge of The Residue2 :[-1] ").strip()
        m2 = input("Multiplicty of The Residue2 :[1] ").strip()
        if not q1:
            q1 = 1
        if not m1:
            m1 = 1
        if not q2:
            q2 = -1
        if not m2:
            m2 = 1
        return q1,m1,q2,m2

    def save(self,s,d,deltax,deltay,deltaz):
        d_name = d*10
        res_name = self.res_name
        atom_name = self.atom_name
        atom_n = self.atom_n
        res_number = self.res_number
        cell = self.cell
        q1,m1,q2,m2 = self.input_info
        if not os.path.exists(f'{d_name:.1f}'):
            os.makedirs(f'{d_name:.1f}')
        with open(f'{d_name:.1f}/{d_name:.1f}.gro', 'w') as f_gro, open(f'{d_name:.1f}/{d_name:.1f}.psi', 'w') as f_psi:
            f_gro.write("%s" %self.blank)
            f_gro.write("%d\n" %self.atoms)
            f_psi.write('molecule sapt_calc {\n')
            f_psi.write(f'{q1} {m1}\n')
            for idx, item in enumerate(zip(self.x,self.y,self.z)) :
                line = f'{res_name[idx]}{atom_name[idx]:>5s}{int(atom_n[idx]):>5d}'
                if idx < res_number[0]:
                    X = item[0] - deltax*s/2
                    Y = item[1] - deltay*s/2
                    Z = item[2] - deltaz*s/2
                    name = ''.join([i for i in atom_name[idx] if not i.isdigit()])  #remove numbers from atom_name
                    f_psi.write(f'{name}{X*10:>8.3f}{Y*10:>8.3f}{Z*10:>8.3f}\n')
                    f_gro.write(f'{line}{X:>8.3f}{Y:>8.3f}{Z:>8.3f}\n')
                else:
                    if idx == res_number[0]:
                        f_psi.write('--\n')
                        f_psi.write(f'{q2} {m2}\n')

                    X = item[0] + deltax*s/2
                    Y = item[1] + deltay*s/2
                    Z = item[2] + deltaz*s/2
                    name = ''.join([i for i in atom_name[idx] if not i.isdigit()])  #remove numbers from atom_name
                    f_psi.write(f'{name}{X*10:>8.3f}{Y*10:>8.3f}{Z*10:>8.3f}\n')
                    f_gro.write(f'{line}{X:>8.3f}{Y:>8.3f}{Z:>8.3f}\n')
            f_psi.write("units angstrom\nno_reorient\nsymmetry c1 }\nset basis aug-cc-pvdz\nenergy('sapt0')\n")
            f_gro.write(cell)
        
        #Gromacs Single Point
        os.system('gmx grompp -f e.mdp -c '+f'{d_name:.1f}/{d_name:.1f}.gro'+' -p topol.top -o '+f'{d_name:.1f}/e.tpr > {d_name:.1f}/gromacs.log')
        cwd = os.getcwd()
        os.chdir(cwd + f'/{d_name:.1f}')
        os.system('gmx mdrun -s e.tpr >> gromacs.log')
        #Sapt PSI4 procedute
        os.system(f'psi4 {d_name:.1f}.psi >> psi4.log')
        #data extraction
        TOTAL, Q, LJ, psi_TOTAL, psi_Q, psi_LJ= self.extract('md.log',f'{d_name:.1f}.out')
        os.chdir(cwd)
        self.data.append([ d*10, TOTAL, Q, LJ, psi_TOTAL, psi_Q, psi_LJ]) #kj, kcal
    def save_all(self):
        self.input_info = self.ask_charge()
        d=self.norm
        s=d/self.n
        shift = s
        self.save(s,d,0,0,0) # initial geometry (probably minimum)
        for j in range(int(self.n*0.8)): #positive shift less points are required
            d=d+s #for the name
            self.save(shift,d,self.deltax,self.deltay,self.deltaz)
            shift += s #for the shift
        d=self.norm
        s=d/self.n
        shift = s
        for j in range(self.n-1): #negative shift
            d=d-s #for the name
            self.save(-shift,d,self.deltax,self.deltay,self.deltaz)
            shift += s
    def extract(self, filename_gro,filename_psi):
        with open(filename_psi) as f_psi, open(filename_gro) as f_gro:
            line_found = False
            while True:
                line_gro = f_gro.readline()
                if "====  A V E R A G E S  ====" in line_gro:
                    line_found = True
                if line_found:
                    for i in range(5):
                        blank = f_gro.readline()
                    header,values = [],[]
                    for i in range(3):
                        header, values = new_line(f_gro,header,values)
                    results = dict(zip(header,values))
                    #Extract values summed
                    LJ = print_sum("LJ", results)
                    Q = print_sum("Coul", results)
                    TOTAL = print_sum("Total", results)
                    break
            while True:
                line_psi = f_psi.readline()
                if "Electrostatics sSAPT0" in line_psi:
                    psi_charge = float(line_psi.split()[4].strip())
                elif "Exchange sSAPT0" in line_psi:
                    psi_exchange = float(line_psi.split()[4].strip())
                elif "Induction sSAPT0" in line_psi:
                    psi_induction = float(line_psi.split()[4].strip())
                elif "Dispersion sSAPT0" in line_psi:
                    psi_disp = float(line_psi.split()[4].strip())
                elif "Total sSAPT0" in line_psi:
                    psi_tot = float(line_psi.split()[4].strip())
                    return TOTAL, Q, LJ, psi_tot, psi_induction+psi_charge, psi_exchange+psi_disp
                    break
    def analysis(self):
        df = self.data
        df = df.sort_values(by = ["Distance"],ignore_index=True)
        df.Total_opls += - df.Total_opls.iloc[-1]
        df.LJ_opls += - df.LJ_opls.iloc[-1]
        df.Coulomb_opls += - df.Coulomb_opls.iloc[-1]
        #kj -> kcal
        df.Total_opls = df.Total_opls/4.184
        df.Coulomb_opls = df.Coulomb_opls / 4.184
        df.LJ_opls = df.LJ_opls / 4.184
        #kcal
        df.Total_sapt += - df.Total_sapt.iloc[-1]
        df.LJ_sapt += - df.LJ_sapt.iloc[-1]
        df.Coulomb_sapt += - df.Coulomb_sapt.iloc[-1]
        self.data = df
        df.to_csv('E.txt', index=False)

def new_line(f,header,values):
    line = f.readline().strip('\n')
    header = header + [line[i:i+15].strip() for i in range(0, len(line), 15)]
    line = f.readline().strip('\n')
    values = values + [line[i:i+15].strip() for i in range(0, len(line), 15)]
    values = [ float(i) for i in values]
    return header, values

def print_sum(string,results):
    matches = [ i for i in results.keys() if string in i]
    return sum([results[i] for i in matches])

obj = coord()
print(obj.data)
