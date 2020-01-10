import math
import random
from matplotlib import pyplot as plt
import glob

def read_tsp_node(path):
    # Open input file
    infile = open(path, 'r')
    # Read instance header
    Name = infile.readline().strip().split()[2] # NAME
    Comment = infile.readline().strip().split()[2] # TYPE
    FileType = infile.readline().strip().split()[2] # COMMENT
    Dimension = infile.readline().strip().split()[2]
    EDGE_WEIGHT_TYPE = infile.readline().strip().split()[2]# DIMENSION
    print('\tFile-Name = ',Name,'.tsp')
    print('\tComment : ',Comment)
    print('\tFileType = ',FileType)
    print('\tDimension = ',Dimension)
    print('\tEDGE_WEIGHT_TYPE = ',EDGE_WEIGHT_TYPE)
    infile.readline()
    
    # Read node list
    nodelist = []
    N = int(Dimension)
    for i in range(0, N):
        x,y = infile.readline().strip().split()[-2:]
        nodelist.append((float(x), float(y)))
    
    infile.close()
    return Name,Comment,FileType,Dimension,nodelist


class ACO:
    class Edge:
        def __init__(self, a, b, weight, ph_init):
            self.a = a
            self.b = b
            self.weight = weight
            self.pherom = ph_init

    class Finder_Root:
        def __init__(self, alpha, beta, num_nodes, edges):
            self.alpha = alpha
            self.beta = beta
            self.num_nodes = num_nodes
            self.edges = edges
            self.tour = None
            self.distance = 0.0

        def node_select(self):
            wheel = 0.0
            nodes_remain = [node for node in range(self.num_nodes) if node not in self.tour]
            total_h = 0.0
            for node_remained in nodes_remain:
                total_h += self.edges[self.tour[-1]][node_remained].weight
            for node_remained in nodes_remain:
                wheel += pow(self.edges[self.tour[-1]][node_remained].pherom, self.alpha) * \
                                  pow((total_h / self.edges[self.tour[-1]][node_remained].weight), self.beta)
            random_value = random.uniform(0.0, wheel)
            w_p = 0.0
            for node_remained in nodes_remain:
                w_p += pow(self.edges[self.tour[-1]][node_remained].pherom, self.alpha) * \
                                  pow((total_h / self.edges[self.tour[-1]][node_remained].weight), self.beta)
                if w_p >= random_value:
                    return node_remained

        def tour_select(self):
            self.tour = [random.randint(0, self.num_nodes - 1)]
            while len(self.tour) < self.num_nodes:
                self.tour.append(self.node_select())
            return self.tour

        def distance_count(self):
            self.distance = 0.0
            for i in range(self.num_nodes):
                self.distance += self.edges[self.tour[i]][self.tour[(i + 1) % self.num_nodes]].weight
            return self.distance

    def __init__(self, method='ACS', size=10, weight=1.0, min_scale=0.001, alpha=1.0, beta=3.0,
                 rho=0.1, w_pherom=1.0, ph_init=1.0, steps=100, nodes=None, labels=None):
        self.method = method
        self.size = size
        self.weight = weight
        self.min_scale = min_scale
        self.rho = rho
        self.w_pherom = w_pherom
        self.steps = steps
        self.num_nodes = len(nodes)
        self.nodes = nodes
        if labels is not None:
            self.labels = labels
        else:
            self.labels = range(1, self.num_nodes + 1)
        self.edges = [[None] * self.num_nodes for _ in range(self.num_nodes)]
        for i in range(self.num_nodes):
            for j in range(i + 1, self.num_nodes):
                self.edges[i][j] = self.edges[j][i] = self.Edge(i, j, math.sqrt(
                    pow(self.nodes[i][0] - self.nodes[j][0], 2.0) + pow(self.nodes[i][1] - self.nodes[j][1], 2.0)),
                                                                ph_init)
        self.Finder_Roots = [self.Finder_Root(alpha, beta, self.num_nodes, self.edges) for _ in range(self.size)]
        self.tour_final = None
        self.distance_final = float("inf")

    def pherom_list(self, tour, distance, weight=1.0):
        ph_to_add = self.w_pherom / distance
        for i in range(self.num_nodes):
            self.edges[tour[i]][tour[(i + 1) % self.num_nodes]].pherom += weight * ph_to_add

    def acs(self):
        for step in range(self.steps):
            for Finder_Root in self.Finder_Roots:
                self.pherom_list(Finder_Root.tour_select(), Finder_Root.distance_count())
                if Finder_Root.distance < self.distance_final:
                    self.tour_final = Finder_Root.tour
                    self.distance_final = Finder_Root.distance
            for i in range(self.num_nodes):
                for j in range(i + 1, self.num_nodes):
                    self.edges[i][j].pherom *= (1.0 - self.rho)

    def elitist(self):
        for step in range(self.steps):
            for Finder_Root in self.Finder_Roots:
                self.pherom_list(Finder_Root.tour_select(), Finder_Root.distance_count())
                if Finder_Root.distance < self.distance_final:
                    self.tour_final = Finder_Root.tour
                    self.distance_final = Finder_Root.distance
            self.pherom_list(self.tour_final, self.distance_final, weight=self.weight)
            for i in range(self.num_nodes):
                for j in range(i + 1, self.num_nodes):
                    self.edges[i][j].pherom *= (1.0 - self.rho)

    def max_min(self):
        for step in range(self.steps):
            tour_iter = None
            dist_iter = float("inf")
            for Finder_Root in self.Finder_Roots:
                Finder_Root.tour_select()
                if Finder_Root.distance_count() < dist_iter:
                    tour_iter = Finder_Root.tour
                    dist_iter = Finder_Root.distance
            if float(step + 1) / float(self.steps) <= 0.75:
                self.pherom_list(tour_iter, dist_iter)
                max_pherom = self.w_pherom / dist_iter
            else:
                if dist_iter < self.distance_final:
                    self.tour_final = tour_iter
                    self.distance_final = dist_iter
                self.pherom_list(self.tour_final, self.distance_final)
                max_pherom = self.w_pherom / self.distance_final
            min_pherom = max_pherom * self.min_scale
            for i in range(self.num_nodes):
                for j in range(i + 1, self.num_nodes):
                    self.edges[i][j].pherom *= (1.0 - self.rho)
                    if self.edges[i][j].pherom > max_pherom:
                        self.edges[i][j].pherom = max_pherom
                    elif self.edges[i][j].pherom < min_pherom:
                        self.edges[i][j].pherom = min_pherom

    def run(self):
        print('Method : {0}'.format(self.method))
        if self.method == 'ACS':
            self.acs()
        elif self.method == 'Elitist':
            self.elitist()
        else: 
            self.max_min()
        print('Root_Result : <- {0} ->'.format(' - '.join(str(self.labels[i]) for i in self.tour_final)))
        print('Distance : {0}\n'.format(round(self.distance_final, 2)))
        return round(self.distance_final, 2)
    
    def plot(self, line_width=1, radius=math.sqrt(2.0), size_ann=8, dpi=120, save=True, name=None,dis=None):
        x = [self.nodes[i][0] for i in self.tour_final]  
        x.append(x[0])
        y = [self.nodes[i][1] for i in self.tour_final]
        y.append(y[0]) 
        plt.plot(x, y, linewidth=line_width, color='red')
        plt.scatter(x, y, s=math.pi * (radius ** 2.0))
        #plt.title(name + '-' + self.method)
        plt.title('ACO Variant: ' + self.method + '      ' + 'TSP Name: ' + name)
        plt.text(0.8,
         0.05, 
         "Distance = {}".format(dis), 
         transform=plt.gca().transAxes) 
        for i in self.tour_final:
            plt.annotate(self.labels[i], self.nodes[i], size=size_ann)
        if save: 
            #if name is None:Should i copy the path here ? 
            name = '/home/matthewvella98/Developer/Optimisation/ACO-TSP/Result-Plot/' + name + '-{0}.png'.format(self.method)
            plt.savefig(name, dpi=dpi) 
        plt.show()
        plt.gcf().clear()



if __name__ == '__main__':
    size = 5
    steps = 50
    in_path = "/home/matthewvella98/Developer/Optimisation/ACO-TSP/Tsp-Data/"
    k = 1
    for path in glob.glob(in_path + "*.tsp"):
        print(" ** Step : ",k,' ** ')
        Name,Comment,FileType,Dimension,node = read_tsp_node(path)
        #print(_nodes)
        Acs = ACO(method='ACS', size=size, steps=steps, nodes=node)
        distance = Acs.run()
        Acs.plot(name = Name,dis = distance )
        Elitist = ACO(method='Elitist',size=size, steps=steps, nodes=node)
        distance = Elitist.run()
        Elitist.plot(name = Name,dis = distance)
        Max_Min = ACO(method='MaxMin', size=size, steps=steps, nodes=node)
        distance = Max_Min.run()
        Max_Min.plot(name = Name,dis = distance)
        k += 1