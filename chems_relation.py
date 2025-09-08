import csv 
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations
from pyvis.network import Network

MAX_DATASET = 1000

def main():
    chem_data = read_csv(MAX_DATASET)

    # create objects
    chem_objs = []
    for chem_name in chem_data.keys():
        elements = chem_data[chem_name]
        chem_obj = chemical(chem_name, elements)
        chem_objs.append(chem_obj)

    #build_network_graph(chem_objs)
    #build_net_graph_v2(chem_objs)
    #build_net_graph_v3(chem_objs)
    build_net_graph_v4(chem_objs)
    print()



def parse_formula(formula):
    # RDKit can parse formulas from molecule objects, but simpler:
    import re
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    elements = []
    for elem, count in tokens:
        elements.extend([elem] * (int(count) if count else 1))
    return sorted(set(elements))

def read_csv(MAX_DATASET):
    
    data = {}
    counter = 0

    filepath = r"C:\Users\adamp\Downloads\unique-reduced-formula.csv"
    
    with open(filepath, 'r', encoding="utf-8") as csvfile:
        csv_reader =  csv.reader(csvfile)
        
        next(csv_reader) # skip headers # ['calc_id', 'reduced_formula', 'count']
        
        for row in csv_reader:
            if counter < MAX_DATASET:
                chem_compound = row[1]
                elements = parse_formula(row[1]) # returns set

                if not chem_compound in data.keys():
                    data[chem_compound] = {}
                
                data[chem_compound] = elements
                counter += 1
            else:
                break

    return data

def build_network_graph(chem_objs):
    
    """networkx 2D graph 
    """
    # Build graph
    G = nx.Graph()

    for chem in chem_objs:
        # Add compound node
        G.add_node(chem.name, type="compound")
        for elem in chem.elements:
            # Add element node
            G.add_node(elem, type="element")
            # Connect compound to element
            G.add_edge(chem.name, elem)

    # Calculate node sizes based on degree
    degree_dict = dict(G.degree())
    node_sizes = [degree_dict[node] * 500 for node in G.nodes()]

    # Color compounds differently from elements
    node_colors = ["skyblue" if G.nodes[node]["type"] == "compound" else "lightgreen" 
                for node in G.nodes()]

    # Draw graph
    plt.figure(figsize=(10, 7))
    pos = nx.spring_layout(G, seed=42)  # nice layout
    nx.draw_networkx(
        G,
        pos,
        with_labels=True,
        node_size=node_sizes,
        node_color=node_colors,
        font_size=10,
        edge_color="gray",
        alpha=0.8
    )
    plt.title("Chemical Relationships", fontsize=14)
    plt.axis("off")
    plt.show()

    print()

def build_net_graph_v2(chem_objs):
    """networkx 2D graph without chem coumpound names
    """
    # Build graph
    G = nx.Graph()

    for chem in chem_objs:
        # Add edges between every pair of elements in the compound
        for e1, e2 in combinations(chem.elements, 2):
            if G.has_edge(e1, e2):
                # Strengthen connection if already exists
                G[e1][e2]["weight"] += 1
            else:
                G.add_edge(e1, e2, weight=1)

    # Calculate node sizes based on degree
    degree_dict = dict(G.degree())
    node_sizes = [degree_dict[node] * 1000 for node in G.nodes()]

    # Draw graph
    plt.figure(figsize=(8, 6))
    pos = nx.spring_layout(G, seed=42)

    # Edge thickness based on weight (# of shared compounds)
    weights = [G[u][v]['weight'] for u, v in G.edges()]

    nx.draw_networkx(
        G,
        pos,
        with_labels=True,
        node_size=node_sizes,
        node_color="lightgreen",
        font_size=10,
        width=weights,
        edge_color="gray",
        alpha=0.8
    )

    plt.title("Element Co-occurrence in Compounds", fontsize=14)
    plt.axis("off")
    plt.show()
    print()

def build_net_graph_v3(chem_objs):
    """Pyvis graph with interactive nodes
    """
    # Build graph
    G = nx.Graph()

    for chem in chem_objs:
        for e1, e2 in combinations(chem.elements, 2):
            if G.has_edge(e1, e2):
                G[e1][e2]["weight"] += 1
            else:
                G.add_edge(e1, e2, weight=1)

    # Convert to Pyvis network
    net = Network(height="700px", width="100%", bgcolor="#222222", font_color="white")

    degree_dict = dict(G.degree())

    for node in G.nodes():
        net.add_node(
            node,
            # Trick: use icon mode to draw text inside node
            shape="circle",
            title=node,   # tooltip
            label="",
            size=degree_dict[node] * 25,
            color="lightgreen",
            font={"size": 0},   # hide default label
            icon={"face": "monospace", "code": node, "size": 50, "color": "black"}
        )

    for u, v, data in G.edges(data=True):
        net.add_edge(u, v, value=data["weight"])

    net.barnes_hut()

    net.write_html("chemicals_network.html")
    print()

def build_net_graph_v4(chem_objs):
    """Pyvis graph with pop up 
    """
    
    # Build graph
    G = nx.Graph()

    for chem in chem_objs:
        for e1, e2 in combinations(chem.elements, 2):
            if G.has_edge(e1, e2):
                G[e1][e2]["weight"] += 1
            else:
                G.add_edge(e1, e2, weight=1)

    # Convert to Pyvis network
    net = Network(height="700px", width="100%", bgcolor="#222222", font_color="white")

    degree_dict = dict(G.degree())

    for node in G.nodes():
        relationships = degree_dict[node]
        # Scaled size: base + multiplier
        size = 15 + relationships * 5
        net.add_node(
            node,
            label=node,
            size=size,
            color="lightgreen",
            title=f"<b>{node}</b><br>Relationships: {relationships}"
        )

    for u, v, data in G.edges(data=True):
        net.add_edge(u, v, value=data["weight"])

    net.barnes_hut()

    net.write_html("chemicals_network.html")

class chemical:
    # Constructor method (called when creating a new object)
    def __init__(self, chem_name, elements):
        self.name = chem_name   # instance attribute
        self.elements = elements     # instance attribute

if __name__ == "__main__":
    main()