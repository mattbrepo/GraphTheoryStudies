import warnings
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import sys
import DetourMatrix

class geeks: 
  def __init__(self, name, roll): 
    self.name = name 
    self.roll = roll

def get_graph(name, connections):
  # extract nodes from graph
  nodes = set([n1 for n1, n2 in connections] + [n2 for n1, n2 in connections])

  # create networkx graph
  G = nx.Graph()

  # add nodes
  for node in nodes:
    G.add_node(node)

  # add edges
  for edge in connections:
    G.add_edge(edge[0], edge[1])

  G.name = name
  return G

def draw_graph(G, file_path):
  f = plt.figure()
  px = 1 / plt.rcParams['figure.dpi']
  plt.subplots(figsize=(300*px, 300*px))
  plt.axis('off')
  plt.tight_layout()
  
  # draw graph
  pos = nx.shell_layout(G)
  nx.draw(G, pos, with_labels=True, node_color='#00b4d9')
  plt.savefig(file_path)
  plt.close('all')

def scipy_to_numpy(scipy_matrix):
  np_matrix = np.array(scipy_matrix.toarray())
  return np_matrix

def draw_matrix(np_matrix, file_path, width=200):
  f = plt.figure()
  px = 1 / plt.rcParams['figure.dpi']
  plt.subplots(figsize=(width*px, 150*px))
  plt.table(cellText=np_matrix, cellLoc='center', loc='center')
  plt.axis('off')
  plt.tight_layout()
  #plt.subplots_adjust(left=0.0, right=1, top=1, bottom=0.0)
  plt.savefig(file_path) #, bbox_inches='tight')
  plt.close('all')

def write_html_matrix(g, M, M_name, f, width=200):
  file_path = 'graphs_images/' + M_name + '_' + g.name + '.png'
  draw_matrix(M, file_path, width)
  f.write('<td><img src="' + file_path + '"></td>')

def getEigens(M, tab_round):
  eigM, eigvectM = np.linalg.eig(M)
  eigM = np.transpose(np.asmatrix(np.around(eigM.real, tab_round)))
  eigvectM = np.around(eigvectM.real, tab_round)
  return eigM, eigvectM

def create_output(gs, GRef):
  f = open("graphs.html", "w")
  f.write('<html><head><style>table, th, td { border: 1px solid black; border-collapse: collapse; }</style></head><body><table><tr>')
  f.write('<th>Name</th>')
  f.write('<th>Plot</th>')
  f.write('<th>Adjacency</th>')
  f.write('<th>eig A</th>')
  f.write('<th>eig vect A</th>')
  f.write('<th>Laplacian</th>')
  f.write('<th>eig L</th>')
  f.write('<th>eig vect L</th>')
  f.write('<th>Normalized Laplacian</th>')
  f.write('<th>eig NL</th>')
  f.write('<th>eig vect NL</th>')
  f.write('<th>Bethe Hessian</th>')
  f.write('<th>eig BH</th>')
  f.write('<th>eig vect BH</th>')
  f.write('<th>Incidence</th>')
  f.write('<th>Detour</th>')
  f.write('<th>eig Dt</th>')
  f.write('<th>eig vect Dt</th>')
  f.write('<th>Detour Ref</th>')
  f.write('<th>eig DtR</th>')
  f.write('<th>eig vect DtR</th>')
  f.write('</tr>')

  eig_table_widh = 100
  tab_round = 2
  descDetMatRef = DetourMatrix.CalcDetour(Gref)
  DtRef = descDetMatRef().astype(int)
  for g in gs:
    f.write('<tr>')
    f.write('<td>' + g.name + '</td>')
    A = scipy_to_numpy(nx.adjacency_matrix(g))
    descDetMat = DetourMatrix.CalcDetour(g)
    Dt = descDetMat().astype(int)
    B = scipy_to_numpy(nx.incidence_matrix(g)).astype(int)
    L = scipy_to_numpy(nx.laplacian_matrix(g))
    NL = scipy_to_numpy(nx.normalized_laplacian_matrix(g))
    NL = np.around(NL, 3)
    BH = np.around(scipy_to_numpy(nx.bethe_hessian_matrix(g)), tab_round)
    
    # eigA = np.transpose(np.asmatrix(np.around(nx.adjacency_spectrum(g).real, tab_round)))
    eigA, eigvectA = getEigens(A, tab_round)
    eigDt, eigvectDt = getEigens(Dt, tab_round)

    #eigL = np.transpose(np.asmatrix(np.around(nx.laplacian_spectrum(g).real, tab_round)))
    #eigNL = np.transpose(np.asmatrix(np.around(nx.normalized_laplacian_spectrum(g).real, tab_round)))
    #eigBH = np.transpose(np.asmatrix(np.around(nx.bethe_hessian_spectrum(g).real, tab_round)))
    eigL, eigvectL = getEigens(L, tab_round)
    eigNL, eigvectNL = getEigens(NL, tab_round)
    eigBH, eigvectBH = getEigens(BH, tab_round)
    
    DtR = DtRef - Dt
    eigDtR, eigvectDtR = getEigens(DtR, tab_round)

    file_path = 'graphs_images/G_' + g.name + '.png'
    draw_graph(g, file_path)
    f.write('<td><img src="' + file_path + '"></td>')

    write_html_matrix(g, A, 'A', f)
    write_html_matrix(g, eigA, 'eigA', f, eig_table_widh)
    write_html_matrix(g, eigvectA, 'eigvectA', f)

    write_html_matrix(g, L, 'L', f)
    write_html_matrix(g, eigL, 'eigL', f, eig_table_widh)
    write_html_matrix(g, eigvectL, 'eigvectL', f)

    write_html_matrix(g, NL, 'NL', f)
    write_html_matrix(g, eigNL, 'eigNL', f, eig_table_widh)
    write_html_matrix(g, eigvectNL, 'eigvectNL', f)

    write_html_matrix(g, BH, 'BH', f)
    write_html_matrix(g, eigBH, 'eigBH', f, eig_table_widh)
    write_html_matrix(g, eigvectBH, 'eigvectBH', f)

    write_html_matrix(g, B, 'B', f)
    #write_html_matrix(g, eigB, 'eigB', f)
    
    write_html_matrix(g, Dt, 'Dt', f)
    write_html_matrix(g, eigDt, 'eigDt', f, eig_table_widh)
    write_html_matrix(g, eigvectDt, 'eigvectDt', f)

    write_html_matrix(g, DtR, 'DtR', f)
    write_html_matrix(g, eigDtR, 'eigDtR', f, eig_table_widh)
    write_html_matrix(g, eigvectDtR, 'eigvectDtR', f)

    f.write('</tr>')
  f.write('</tr></table></body></html>')
  f.close()

#
# Main
#

warnings.simplefilter(action='ignore', category=FutureWarning)

if False:
  #connections = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]
  G = get_graph("1", [(0, 3), (1, 3), (2, 3)])
  print(G.name)
  A = scipy_to_numpy(nx.adjacency_matrix(G))
  descDetMat = DetourMatrix.CalcDetour(G) # DetourMatrix.DetourMatrixCache()
  Dt = descDetMat()
  eigA = nx.adjacency_spectrum(G).real
  B = scipy_to_numpy(nx.incidence_matrix(G))
  #eigB = np.transpose(np.asmatrix(np.around(np.linalg.eigvals(B).real, 3)))

  eigDt, eigvectDt = np.linalg.eig(Dt)
  print(eigvectDt)
  #draw_graph(G, 'G.png')
  #draw_matrix(A, 'A.png')
  #draw_matrix(Dt, 'Dt.png')
  print(B)
  sys.exit()

# 3 vertices
if False:
  list = []
  list.append(get_graph("1", [(0, 2), (1, 2)]))
  list.append(get_graph("2", [(0, 1), (0, 2), (1, 2)]))
  
  Gref = get_graph("ref", [(0, 1), (0, 2), (1, 2)])
  create_output(list, Gref)
  sys.exit()

# 4 vertices
if True:
  list = []
  list.append(get_graph("1", [(0, 3), (1, 3), (2, 3)]))
  list.append(get_graph("2", [(0, 2), (0, 3), (1, 3)]))
  list.append(get_graph("3", [(0, 2), (0, 3), (1, 3), (2, 3)]))
  list.append(get_graph("4", [(0, 2), (0, 3), (1, 2), (1, 3)]))
  list.append(get_graph("5", [(0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]))
  list.append(get_graph("6", [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]))
  
  Gref = get_graph("ref", [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])
  create_output(list, Gref)
  sys.exit()
