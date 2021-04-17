import re
from anytree import Node, RenderTree, PreOrderIter
from anytree.exporter import  DotExporter
import copy
import os


def find_by_order(matches_numbering):

    # A very high number in place of infinity to use as an initial value
    lowest = 9999999999
    lowest_index = 0

    for x in matches_numbering:
        if int(x) < lowest:
            lowest = int(x)
            lowest_index = matches_numbering.index(x)
            # print(lowest_index)

    return lowest_index


def wrong_material(x, text_only):
    """
    Manually checks chemicals that are of no interest to us.

    :param x: the starting material entry
    :param text_only: the text file containing the unannotated text
    :return: false if the entry is not simply called "Compound X", true if the entry is called this, making it useless.
    """

    matches_wrong1 = re.search("compound\s+" + re.escape(x) +"", text_only)
    matches_wrong2 = re.search("Compound\s+" + re.escape(x) +"", text_only)

    if matches_wrong1 is None and matches_wrong2 is None:
        return False
    return True


def edgeattrfunc(node,child):
    return 'label="%s"' % (child.origin)


class LinkFinder:
    def __init__(self, directory, start_tag="STARTING_MATERIAL", end_tag="REACTION_PRODUCT", intermediate_tags=["WORKUP", "REACTION_STEP"], start_number=0, end_number=10000, number_length=4, debug=False):
        """
        Class for finding and searching through cascading reaction chains
        :param directory: Directory to find the files in
        :param start_tag: Entity tag to use as start of reaction
        :param end_tag: Entity tag to use as end of reaction
        :param intermediate_tags: Tags to use as intermediates
        :param start_number: Number of first file
        :param end_number: Number of last file
        :param number_length: Character length of the file numbers
        :param debug: Turn debug print statements on or off
        """
        self.directory = directory
        self.startTag = start_tag
        self.endTag = end_tag
        self.intermediateTags=intermediate_tags
        self.startNumber = start_number
        self.endNumber = end_number
        self.numberLength = number_length
        self.debug = debug

        self.foundLinks = {}
        self.chain = {}
        self.root = None
        self.rootTruncated = None

    def __relation_link(self, forward, ending, covered, matches_intermediate, text):

        """
        Finds whether two entities are connected by a series of relations in any way through intermediates

        :param forward: the entity being checked
        :param ending: the target entity
        :param covered: a list of entities already checked
        :param matches_intermediate: a list of entities that could potentially act as intermediates
        :param text: the annotation text being read from
        :return: True if there is a link found, False if there is not
        """

        # print("finding connection between "+forward+" and " +ending)

        relationEnd = re.findall("R\d+(?=\tARG\w\sArg1:" + forward + "\sArg2:" + ending + "\s)", text)
        if relationEnd:
            # print("connection found between "+forward+" and "+ending)
            return True
        for step in matches_intermediate:
            relation1 = re.findall("R\d+(?=\tARG\w\sArg1:" + step + "\sArg2:" + forward + "\s)", text)
            relation2 = re.findall("R\d+(?=\tARG\w\sArg1:" + forward + "\sArg2:" + step + "\s)", text)
            if relation1 and not relation1[0] in covered:
                covered.append(relation1[0])
                if self.__relation_link(step, ending, covered, matches_intermediate, text):
                    return True
            if relation2 and not relation2[0] in covered:
                covered.append(relation2[0])
                if self.__relation_link(step, ending, covered, matches_intermediate, text):
                    return True
        # print("failed to find connection between "+forward+" and "+ ending)
        return False

    def __recursive_tree_form(self, forward, root):

        """
        Creates a tree onto the node based on the known start and end connections

        :param forward: The node onto which links will be attached
        :param root: The root node
        :return: Always returns True, modifies the inputted nodes by attaching children
        """

        if forward.name in self.chain and (self.chain[forward.name].parent == forward.parent or self.chain[forward.name].parent == root or forward.parent == root):
            # print("in chain, depth: " + str(forward.depth) + " vs. " + str(chain[forward.name].depth))
            if self.chain[forward.name].depth < forward.depth:
                # print("toast")
                forward.children = self.chain[forward.name].children
                self.chain[forward.name].parent = None
                self.chain[forward.name] = forward
            else:
                forward.parent = None
        else:
            self.chain[forward.name] = forward
            for entry in self.foundLinks:
                value = self.foundLinks[entry]
                for match in value:
                    if match[0] == forward.name:
                        newNode = Node(match[1], parent=forward, origin=entry)
                        self.__recursive_tree_form(newNode, root)

        return True

    def link_files(self):

        """
        Runs the link finding system through the files and arranges them into a tree structure
        :return: True if the process worked successfully, False if no links were found at all
        """

        file_number = 0
        end = False
        while not end:
            try:
                # Open both the annotation file and the text file
                file_ann = open(self.directory + str(file_number).zfill(4) + ".ann", "r")
                file_text = open(self.directory + str(file_number).zfill(4) + ".txt","r")

                if self.debug:
                    print(file_ann.name)

                ann_text = file_ann.read()
                content_text = file_text.read()

                # Finds all entities that match the starting tag
                matches_start = re.findall("T\d+(?=\t"+self.startTag+")", ann_text)
                # Finds the number of the entity in the order
                matches_start_numbering = re.findall("(?<=\t"+self.startTag+"\s)\d+",ann_text)

                # Finds all entities that match the ending tag
                matches_end = re.findall("T\d+(?=\t"+self.endTag+")", ann_text)
                # Finds the number of the entity in the order
                matches_end_numbering = re.findall("(?<=\t"+self.endTag+"\s)\d+",ann_text)

                # Finds all entities that match the intermediate tags
                matches_intermediate = []
                for intermediate in self.intermediateTags:
                    matches = re.findall("T\d+(?=\t"+intermediate+")", ann_text)
                    for match in matches:
                        matches_intermediate.append(match)

                # Initializes the earliest occurence of and end product to a dummy value
                earliest_end = "not correct"

                # Sets the earliest end and start entities to their correct values
                if matches_start and matches_end:
                    earliest_start = matches_start[find_by_order(matches_start_numbering)]
                    earliest_end = matches_end[find_by_order(matches_end_numbering)]

                # Iterates through start and end tags to find which are connected, and adds those links to foundLinks
                for x in matches_start:
                    for y in matches_end:

                        # Checks if the start and end are related, or if they are the first two in the file
                        if (self.__relation_link(x, y, [], matches_intermediate, ann_text) or (x == earliest_start and y == earliest_end)):
                            if self.debug:
                                print(x)
                                print(y)

                            # Gets the textual name of the entity
                            starting_tag = re.findall("" + x + "\t"+self.startTag+"\s\d+\s\d+\t.+", ann_text)
                            starting_tag = starting_tag[0]
                            starting_tag = re.sub("" + x + "\t"+self.startTag+"\s\d+\s\d+\t", "", starting_tag)

                            ending_tag = re.findall("" + y + "\t"+self.endTag+"\s\d+\s\d+\t.+", ann_text)
                            ending_tag=ending_tag[0]
                            ending_tag = re.sub("" + y + "\t"+self.endTag+"\s\d+\s\d+\t", "", ending_tag)

                            # Checks that the entities contain useful information
                            if not wrong_material(starting_tag, content_text):

                                # Adds the given start and end pair to foundLinks
                                material_io = [starting_tag, ending_tag]
                                try:
                                    self.foundLinks[file_number].append(material_io)
                                except KeyError:
                                    self.foundLinks[file_number] = []
                                    self.foundLinks[file_number].append(material_io)

                file_ann.close()
                file_number = file_number + 1

            except FileNotFoundError:
                # Checks if the attempted file is over the maximum file number
                if file_number > self.endNumber:

                    # Checks if any links were found at all
                    if self.foundLinks:
                        if self.debug:
                            print(self.foundLinks)
                        end = True
                    else:
                        print("Failed to find a match!")
                        return False
                # Simply increments the file number to skip over gaps in the number sequence
                file_number = file_number + 1

        root = Node("root",origin="-1")

        for entry in self.foundLinks:
            value = self.foundLinks[entry]
            for match in value:

                nodule = Node(match[0], parent=root, origin=entry)
                if self.debug:
                    print(match)
                    print(nodule.name)
                    print(nodule)
                self.__recursive_tree_form(nodule, root)
        self.root = root

        root_truncated = copy.deepcopy(root)
        ref_node = Node("Success", parent=root_truncated, origin=-1)
        for node in PreOrderIter(root_truncated):
            if node.depth == 1:
                ref_node.parent = None
                ref_node = node

            if node.depth > 2:
                if self.debug:
                    print(node)
                    print("is level 3")
                ref_node = Node("Success", parent=root_truncated, origin=-1)

        self.rootTruncated = root_truncated

        return True

    def ascii_export(self, truncated=True):
        if truncated:
            print(RenderTree(self.rootTruncated))
        else:
            print(RenderTree(self.root))

    def export_svg(self, file_name, truncated=True):

        if truncated:
            DotExporter(self.rootTruncated, edgeattrfunc=edgeattrfunc).to_dotfile(file_name + ".dot")
        else:
            DotExporter(self.root, edgeattrfunc=edgeattrfunc).to_dotfile(file_name + ".dot")
        old_file = open(file_name + ".dot", "r")
        new_file = open(file_name + ".new.dot", "w")

        new_lines = []
        for line in old_file.readlines():
            if '"root"' not in line:
                new_lines.append(line)
        new_file.writelines(new_lines)
        new_file.close()
        os.system("dot " + file_name + ".new.dot -Nshape=box -GoutputMode=nodesFirst -Grankdir=LR -T svg -O")

    def find_by_key(self,key):

        for node in PreOrderIter(self.root):
            if key in node.name:
                print(RenderTree(node))

    def find_by_file(self,file_number):
        for node in PreOrderIter(self.root):
            try:
                if node.origin == file_number:
                    print(RenderTree(node))
            except AttributeError:
                print(node)


