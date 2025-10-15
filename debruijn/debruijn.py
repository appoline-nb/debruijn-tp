#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""
'''  -i fichier fastq single end
-k longueur des kmer (optionnel - défaut 22)
 -o fichier output avec les contigs'''


import argparse
import os
import sys
from pathlib import Path
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
    random_layout,
    draw,
    spring_layout,
)
import matplotlib
from operator import itemgetter
import random

random.seed(9001)
from random import randint
import statistics
import textwrap
import matplotlib.pyplot as plt
from typing import Iterator, Dict, List

matplotlib.use("Agg")

__author__ = "Appoline NABI"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Appoline NABI"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "appoline.nabi@gmail.com"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage="{0} -h".format(sys.argv[0])
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()

#Question 1 : 
def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.
    : param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with open(fastq_file) as fasq:
        while True: 
            id_line = fasq.readline().strip() 
            if not id_line: 
                break 
            read = fasq.readline().strip() 
            plus = fasq.readline() 
            quality = fasq.readline() 
            yield read 
    pass

def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.
    : param read: (str) Sequence of a read.
    : return: A generator object that provides the kmers (str) of size kmer_size."""
    for k in range(len(read) - kmer_size + 1):
        yield read[k : k + kmer_size]
    pass


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file
     :param fastq_file: (str) Path to the fastq file.
    : return: A dictionnary object that identify all kmer occurrences.
    """
from pathlib import Path
from typing import Dict

def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]: 
    kmer_dict: Dict[str, int] = {}
    with open(fastq_file) as f:
        while True:
            id_line = f.readline().strip()
            if not id_line:
                break
            read = f.readline().strip()
            f.readline()  
            f.readline()  

            for i in range(len(read) - kmer_size + 1):
                kmer = read[i:i+kmer_size]
                kmer_dict[kmer] = kmer_dict.get(kmer, 0) + 1
    return kmer_dict

    pass

#Question 1 : 
from networkx import DiGraph   
import networkx as nx
def build_graph(kmer_dict: Dict[str, int]) -> DiGraph: 
    """Build the debruijn graph
    :param kmer_dict: A dictionnary object that identify all kmer occurrences.
    :return: A directed graph (nx) of all kmer substring and weight (occurrence).
    """
    G = DiGraph()
    for kmer, count in kmer_dict.items():
        prefix = kmer[:-1]
        suffix = kmer[1:]
        if G.has_edge(prefix, suffix):
            G[prefix][suffix]["weight"] += count
        else:
            G.add_edge(prefix, suffix, weight=count)
    return G 
    pass

#Question 3a:
def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    for path in path_list:
        if not path:
            continue
        
        start_idx = 0 if delete_entry_node else 1
        end_idx   = len(path) if delete_sink_node else len(path) - 1
        
        if end_idx <= start_idx:
            continue
        nodes_to_remove = path[start_idx:end_idx]

        G.remove_nodes_from([n for n in nodes_to_remove if G.has_node(n)])
    return G
    pass

#Question 3b: 
from typing import List
import random
import statistics
import networkx as nx
def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    if not path_list:
        return G
    
    candidates = list(range(len(path_list)))
    if len(weight_avg_list) >= 2 and statistics.stdev(weight_avg_list) > 0:
        max_w = max(weight_avg_list)
        candidates = [i for i, w in enumerate(weight_avg_list) if w == max_w]
    elif len(path_length) >= 2 and statistics.stdev(path_length) > 0:
        max_l = max(path_length)
        candidates = [i for i, l in enumerate(path_length) if l == max_l]
    else:
        random.seed(9001)
        candidates = [random.randint(0, len(path_list) - 1)]
    if len(candidates) > 1:
        random.seed(9001)
        keep_idx = random.choice(candidates)
    else:
        keep_idx = candidates[0]
    to_remove = [p for i, p in enumerate(path_list) if i != keep_idx]
    G = remove_paths(G, to_remove, delete_entry_node, delete_sink_node)
    return G
    pass


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )

#Question 3b :
def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    if not nx.has_path(G, ancestor_node, descendant_node):
            return G

    paths = list(nx.all_simple_paths(G, source=ancestor_node, target=descendant_node))
    if len(paths) <= 1:
        return G
    lengths = [len(p) for p in paths]
    weights = []
    for path in paths:
        edge_weights = []
        for u, v in zip(path, path[1:]):
            if G.has_edge(u, v):
                edge_weights.append(G[u][v].get("weight", 1))
        weights.append(sum(edge_weights) / len(edge_weights) if edge_weights else 0.0)

    return select_best_path(
        G,
        path_list=paths,
        path_length=lengths,
        weight_avg_list=weights,
        delete_entry_node=False,
        delete_sink_node=False,
    )
    pass

#Question 3b : 
def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles
    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    changed = True
    while changed:
        changed = False

        for node in list(G.nodes()):
            preds = list(G.predecessors(node))
            if len(preds) < 2:
                continue
            lca_found = None
            for i in range(len(preds)):
                for j in range(i + 1, len(preds)):
                    lca = nx.lowest_common_ancestor(G, preds[i], preds[j])
                    if lca is not None:
                        lca_found = lca
                        break
                if lca_found is not None:
                    break
            if lca_found is None:
                continue
            before = (G.number_of_nodes(), G.number_of_edges())
            G2 = solve_bubble(G, lca_found, node)
            after = (G2.number_of_nodes(), G2.number_of_edges())
            if after != before:
                G = G2
                changed = True
                break  
    return G
    pass

#Question 4 :
def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    changed = True

    while changed:
        changed = False

        current_starts = [s for s in starting_nodes if G.has_node(s)]

        for node in list(G.nodes()):
            if not G.has_node(node):
                continue
            reachable_starts = []
            for s in current_starts:
                if s != node and G.has_node(s) and nx.has_path(G, s, node):
                    reachable_starts.append(s)
            if len(reachable_starts) < 2:
                continue

            paths = []
            for s in reachable_starts:
                paths += list(nx.all_simple_paths(G, source=s, target=node))
            if len(paths) <= 1:
                continue

            lengths = [len(p) for p in paths]
            weights = []
            for p in paths:
                ew = [G[u][v].get("weight", 1) for u, v in zip(p, p[1:]) if G.has_edge(u, v)]
                weights.append(sum(ew) / len(ew) if ew else 0.0)

            G2 = select_best_path(
                G,
                path_list=paths,
                path_length=lengths,
                weight_avg_list=weights,
                delete_entry_node=True,   
                delete_sink_node=False,   
            )

            if (G2.number_of_nodes(), G2.number_of_edges()) != (G.number_of_nodes(), G.number_of_edges()):
                G = G2
                changed = True
                break

    return G
    pass

#Question 4:
def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    G = graph.copy()
    changed = True

    while changed:
        changed = False

        for node in list(G.nodes()):
            if not G.has_node(node):
                continue
            if G.out_degree(node) < 2:   
                continue
            ends_here = [t for t in ending_nodes if t != node and G.has_node(t) and nx.has_path(G, node, t)]
            if len(ends_here) < 2:
                continue
            paths = []
            for t in ends_here:
                paths += list(nx.all_simple_paths(G, source=node, target=t))
            if len(paths) < 2:
                continue
            lengths = [len(p) for p in paths]
            weights = []
            for p in paths:
                ew = [G[u][v].get("weight", 1) for u, v in zip(p, p[1:]) if G.has_edge(u, v)]
                weights.append(sum(ew) / len(ew) if ew else 0.0)
            G2 = select_best_path(
                G,
                path_list=paths,
                path_length=lengths,
                weight_avg_list=weights,
                delete_entry_node=False,   
                delete_sink_node=True,
            )     
            if (G2.number_of_nodes(), G2.number_of_edges()) != (G.number_of_nodes(), G.number_of_edges()):
                G = G2
                changed = True
                break  
    return G
    pass

#Question 2 :
def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    starting_nodes = []
    for node in graph.nodes():
        if len(list(graph.predecessors(node))) == 0:
            starting_nodes.append(node)
    return starting_nodes
    pass

#Question 2:
def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    sink_nodes = []
    for node in graph.nodes():
        if len(list(graph.successors(node))) == 0:
            sink_nodes.append(node)
    return sink_nodes
    pass

#Question 2:
def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    for start in starting_nodes:
        for end in ending_nodes:
            if nx.has_path(graph, start, end):
                for path in nx.all_simple_paths(graph, source=start, target=end):
                    seq = path[0]
                    for node in path[1:]:
                        seq += node[-1]
                    contigs.append((seq, len(seq)))

    return contigs
    pass

#Question 2:
def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contig_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as f:
        for i, (contig, length) in enumerate(contigs_list):
            f.write(f">contig_{i} len={length}\n")
            f.write(textwrap.fill(contig, width=80))
            f.write("\n")
    pass


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]

    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )

    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    args = get_arguments()

    # 1) Lecture + construction du graphe
    print("[1/5] Construction du graphe")
    kmer_dict = build_kmer_dict(args.fastq_file, args.kmer_size)
    graph = build_graph(kmer_dict)
    # dessin du graphe brut
    if args.graphimg_file:
        try:
            draw_graph(graph, args.graphimg_file)
            print(f"    Graphe brut sauvegardé dans: {args.graphimg_file}")
        except Exception as e:
            print(f"    (warning) Impossible de dessiner le graphe: {e}")

    # 2) Résolution des bulles
    print("[2/5] Simplification des bulles…")
    graph = simplify_bubbles(graph)

    # 3) Résolution des pointes d’entrée et de sortie
    print("[3/5] Nettoyage des pointes…")
    starts = get_starting_nodes(graph)
    graph = solve_entry_tips(graph, starts)

    ends = get_sink_nodes(graph)
    graph = solve_out_tips(graph, ends)

    # Recalcule les extrémités après nettoyage
    starts = get_starting_nodes(graph)
    ends = get_sink_nodes(graph)

    # 4) Extraction des contigs
    print("[4/5] Extraction des contigs…")
    contigs = get_contigs(graph, starts, ends)
    if not contigs:
        print("Aucun contig extrait — vérifie -k et l'entrée FASTQ.")
        sys.exit(1)

    # Écriture de tous les contigs
    args.output_file.parent.mkdir(parents=True, exist_ok=True)
    save_contigs(contigs, args.output_file)
    print(f"    Contigs écrits dans: {args.output_file}")

    # Sauvegarde en plus le plus long contig pour BLAST pratique
    longest_seq, longest_len = max(contigs, key=itemgetter(1))
    longest_path = args.output_file.with_name("contig")
    with open(longest_path, "w") as f:
        f.write(f">contig len={longest_len}\n")
        f.write(textwrap.fill(longest_seq, width=80))
        f.write("\n")
    print(f"    Plus long contig sauvegardé dans: {longest_path} (len={longest_len})")

    # 5) Récap : 
    print("[5/5] Fini")
    print(f"    #contigs: {len(contigs)} | max: {longest_len} nt | k={args.kmer_size}")


if __name__ == "__main__":  # pragma: no cover
    main()