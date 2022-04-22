#!/usr/bin/python3
# -*- coding: utf-8 -*-
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, Circle, Ellipse, FancyBboxPatch
import math
import matplotlib.pyplot as plt
import shutil
import optparse
import pcbnew
import numpy as np
import os
import re
import csv
#import sys
#reload(sys)
#sys.setdefaultencoding("utf-8")

ignore_footpr = ["FIDUCIAL", "^TP[-_]", "^TP$", "^SP$",
                 "^SCREW_[S]?[HM][0-9\.]{1,3}[_W]{0,2}$", "^Jumper[1-9]_Triangle"]
ignore_value = ["FIDUCIAL", "DNP"]

bom_table_items_per_page = 50

bom_csv_open = 0
bom_csv_writer = None

def create_board_bom(pcb, boards, bom_table, bom_table_start, bom_table_end):
	global bom_csv_open, bom_csv_writer
	plt.figure(figsize=(5.8, 8.2))

	color_pad = "lightgray"
	color_pad_highlight = "#990000"
	color_pad_pin1 = "#FF0000"
	color_bbox1 = "None"
	color_bbox2 = "#E9AFAF"

	plt.text(0.5, 1.1, "BOM List: " + os.path.splitext(pcb.GetFileName())[0], wrap=False, horizontalalignment='center', verticalalignment='top')

	if boards == None:
		columns = ('Anzahl', 'Footprint', 'Value')
		column_widths = [0.1, 0.45, 0.45]
	else:
		columns = ('Anzahl x %d' % (boards), 'Anzahl', 'Footprint', 'Value')
		column_widths = [0.1, 0.1, 0.4, 0.4]

	if bom_csv_open == 0 and options.csv_bom == True:
		bom_csv_open = open(os.path.splitext(file)[0] + "_bom.csv", 'w')
		bom_csv_writer = csv.writer(bom_csv_open)
		bom_csv_writer.writerow(['total_qty', 'board_qty', 'footprint', 'value', 'refs'])

	cell_text = []
	for i, bom_row in enumerate(bom_table):
		if i >= bom_table_start and i < (bom_table_start+bom_table_items_per_page):
			qty, value, footpr, smt, highlight_refs = bom_row
			if boards == None:
				cell_text.append([str(qty), footpr, value])
			else:
				cell_text.append([str(qty * boards), str(qty), footpr, value])
			if options.csv_bom == True:
				bom_csv_writer.writerow([str(qty * (1 if boards == None else boards)), str(qty), footpr, value, ", ".join(highlight_refs)])

	#the_table = plt.table(cellText=cell_text, rowLabels=rows, rowColours=colors, colLabels=columns, loc='bottom')
	table = plt.table(cellText=cell_text, colLabels=columns, loc='upper left', bbox=[0.0, 0, 1, 1], colWidths=column_widths)
	table.set_fontsize(12)
	#table.scale(1.5, 1.5)

	plt.axis('off')

def create_board_figure(pcb, bom_row, boards, layer=pcbnew.F_Cu):
    qty, value, footpr, smt, highlight_refs = bom_row

    plt.figure(figsize=(5.8, 8.2))
    ax = plt.subplot("111", aspect="equal")

    color_pad = "lightgray"
    color_pad_highlight = "#990000"
    color_pad_pin1 = "#FF0000"
    color_bbox1 = "None"
    color_bbox2 = "#E9AFAF"

    # get board edges (assuming rectangular, axis aligned pcb)
    edge_coords = []
    for d in pcb.GetDrawings():
        if (d.GetLayer() == pcbnew.Edge_Cuts):
            edge_coords.append(d.GetStart())
            edge_coords.append(d.GetEnd())
    edge_coords = np.asarray(edge_coords) * 1e-6
    board_xmin, board_ymin = edge_coords.min(axis=0)
    board_xmax, board_ymax = edge_coords.max(axis=0)

    # add title
    if not boards or boards == 0:
        ax.text(board_xmin + .5 * (board_xmax - board_xmin), board_ymin - 0.5,
            "%dx %s, %s" % (qty, value, footpr), wrap=True,
            horizontalalignment='center', verticalalignment='bottom')
    else:
        ax.text(board_xmin + .5 * (board_xmax - board_xmin), board_ymin - 0.5,
            "(Total: %dx) %dx %s, %s" % (qty * boards, qty, value, footpr), wrap=True,
            horizontalalignment='center', verticalalignment='bottom')

    # add ref list
    ax.text(board_xmin + .5 * (board_xmax - board_xmin), board_ymax + 0.5,
            str(len(highlight_refs)) + " - " + (", ".join(highlight_refs)), wrap=True,
            horizontalalignment='center', verticalalignment='top')

    # draw parts
    for m in pcb.GetFootprints():
        if m.GetLayer() != layer:
            continue
        ref, center = m.GetReference(), np.asarray(m.GetCenter()) * 1e-6

        # do the same checks and replacements as in bom generation
        mvalue = m.GetValue()
        try:
            mfootpr = str(m.GetFPID().GetFootprintName())
        except:
            mfootpr = str(m.GetFPID().GetLibItemName())

        if options.precise != True:
            mfootpr = re.sub(r'^([A-Z]{1,3}[0-9]{4})[LNM]$', r'\1', mfootpr)

        if options.dontcongroup == True:
            if re.match('^[XM][0-9]', ref):
                mvalue = mfootpr

        if ref in highlight_refs and mvalue == value and mfootpr == footpr:
            highlight = True
        else:
            highlight = False

        # bounding box
        mrect = m.GetBoundingBox(False, False)
        mrect_pos = np.asarray(mrect.GetPosition()) * 1e-6
        mrect_size = np.asarray(mrect.GetSize()) * 1e-6
        if layer == pcbnew.B_Cu:
            mrect_pos[0] = (((board_xmax - board_xmin) - (mrect_pos[0] - board_xmin)) + board_xmin) - mrect_size[0]
        rct = Rectangle(mrect_pos, mrect_size[0], mrect_size[1])
        rct.set_color(color_bbox2 if highlight else color_bbox1)
        rct.set_zorder(-1)
        if highlight:
            rct.set_linewidth(.1)
            rct.set_edgecolor(color_pad_highlight)
        ax.add_patch(rct)

        # center marker
        if highlight:
            if layer == pcbnew.B_Cu:
                center[0] = ((board_xmax - board_xmin) - (center[0] - board_xmin)) + board_xmin
            plt.plot(center[0], center[1], ".",
                     markersize=mrect_size.min(), color=color_pad_highlight)

        # plot pads
        for p in m.Pads():
            pos = np.asarray(p.GetPosition()) * 1e-6
            if layer == pcbnew.B_Cu:
                pos[0] = ((board_xmax - board_xmin) - (pos[0] - board_xmin)) + board_xmin
            size = np.asarray(p.GetSize()) * 1e-6 * .9

            is_pin1 = p.GetPadName() == "1" or p.GetPadName() == "A1" or p.GetPadName() == "C"
            shape = p.GetShape()
            offset = p.GetOffset()  # TODO: check offset

            # pad rect
            angle = p.GetOrientation() * 0.1
            cos, sin = np.cos(
                np.pi / 180. * angle), np.sin(np.pi / 180. * angle)
            dpos = np.dot([[cos, -sin], [sin, cos]], -.5 * size)

            if shape == 1:
                rct = Rectangle(pos + dpos, size[0], size[1], angle=angle)
            elif shape == 2:
                rct = Rectangle(pos + dpos, size[0], size[1], angle=angle)
            elif shape == 4:
                rct = Rectangle(pos + dpos, size[0], size[1], angle=angle)
            elif shape == 5:
                rct = Rectangle(pos + dpos, size[0], size[1], angle=angle)
            elif shape == 0:
                rct = Ellipse(pos, size[0], size[1], angle=angle)
            else:
                print("Unsupported pad shape " + str(shape))
                continue
            rct.set_linewidth(0)
            rct.set_color(color_pad_highlight if highlight else color_pad)
            rct.set_zorder(1)
            # highlight pin1
            if highlight and is_pin1:
                rct.set_color(color_pad_pin1)
                rct.set_linewidth(.1)
                rct.set_edgecolor(color_pad_highlight)
            ax.add_patch(rct)

    plt.xlim(board_xmin, board_xmax)
    plt.ylim(board_ymax, board_ymin)

    plt.axis('off')


def natural_sort(l):
    """
    Natural sort for strings containing numbers
    """
    def convert(text): return int(text) if text.isdigit() else text.lower()

    def alphanum_key(key): return [convert(c)
                                   for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)


def generate_bom(pcb, filter_layer=None):
    """
    Generate BOM from pcb layout.
    :param filter_layer: include only parts for given layer
    :return: BOM table (qty, value, footprint, refs)
    """

    # build grouped part list
    part_groups = {}
    for m in pcb.GetFootprints():
        smd_part = 0
        # filter part by layer
        if filter_layer is not None and filter_layer != m.GetLayer():
            continue
        # group part refs by value and footprint
        value = m.GetValue()
        try:
            footpr = str(m.GetFPID().GetFootprintName())
        except:
            footpr = str(m.GetFPID().GetLibItemName())

        if options.precise != True:
            footpr = re.sub(r'^([A-Z]{1,3}[0-9]{4})[LNM]$', r'\1', footpr)

        # check for smd pads
        for p in m.Pads():
            if p.GetAttribute() == 1:
                smd_part = 1

        reference = m.GetReference()
        if options.dontcongroup == True:
            if re.match('^[XM][0-9]', reference):
                value = footpr

        group_key = (value, footpr, smd_part)
        refs = part_groups.setdefault(group_key, [])
        refs.append(reference)

    # build bom table, sort refs
    bom_table_smd = []
    bom_table_thd = []
    for (value, footpr, smd_part), refs in part_groups.items():
        add_footprint = 1
        for regex in ignore_value:
            m = re.match(regex, value)
            if m:
                add_footprint = 0

        for regex in ignore_footpr:
            m = re.match(regex, footpr)
            if m:
                add_footprint = 0

        line = (len(refs), value, footpr, smd_part, natural_sort(refs))
        if add_footprint == 1 and smd_part != 1:
            bom_table_thd.append(line)
        elif add_footprint == 1:
            bom_table_smd.append(line)
        else:
            print("Ignoring footprint " + footpr + " with value " + value)

    # sort table by reference prefix and quantity
    def sort_func(row):
        qty, _, _, _, rf = row
        ref_ord = {"R": 3, "C": 3, "L": 1, "D": 1,
                   "J": -1, "P": -1}.get(rf[0][0], 0)
        return -ref_ord, -qty

    bom_table_smd = sorted(bom_table_smd, key=sort_func)
    bom_table_thd = sorted(bom_table_thd, key=sort_func)

    bom_table = bom_table_smd + bom_table_thd

    return bom_table


if __name__ == "__main__":
    parser = optparse.OptionParser(
        description='KiCad PCB pick and place assistant')
    parser.add_option('-s', '--split', action="store_true",
                      dest="split", help="split into one file per component")
    parser.add_option('-p', '--precise', action="store_true", dest="precise",
                      help="be precise with footprints (else C0402N == C0402L)")
    parser.add_option('-m', '--dontconnectorgroup', action="store_true",
                      dest="dontcongroup", help="ignore connector values when grouping parts")
    parser.add_option('-f', '--folders', action="store_true", dest="fsort",
                      help="sort pdfs (when splitting) to folders tht smt smt/passives")
    parser.add_option('-n', '--panels', action="store", dest="boards", type="int",
                      help="set number of boards/panels to populate")
    parser.add_option('-b', '--bom-only', action="store_true", dest="bom_only",
                      help="print only bom, then exit")
    parser.add_option('-c', '--csv-bom', action="store_true", dest="csv_bom",
                      help="create csv bom file")
    options, args = parser.parse_args()

    if len(args) != 1:
        parser.error("wrong number of arguments")

    # build BOM
    file = args[0]
    print("Loading %s" % file)
    pcb = pcbnew.LoadBoard(file)
    bom_table = generate_bom(pcb, filter_layer=None)
    bom_table_bot = generate_bom(pcb, filter_layer=pcbnew.B_Cu)
    bom_table_top = generate_bom(pcb, filter_layer=pcbnew.F_Cu)

    # for each part group, print page to PDF
    if options.split != True:
        fname_out = os.path.splitext(file)[0] + "_picknplace.pdf"
        with PdfPages(fname_out) as pdf:
            bom_table_items = len(bom_table)
            if bom_table_items <= bom_table_items_per_page:
                bom_table_pages = 1
            else:
                bom_table_pages = int(math.ceil(bom_table_items / bom_table_items_per_page))
                bom_table_items_per_page = int(math.ceil(bom_table_items / bom_table_pages))
                
            for st in range(0, bom_table_pages):
                create_board_bom(pcb, options.boards, bom_table, st * bom_table_items_per_page, (st+1) * bom_table_items_per_page)
                pdf.savefig()
                plt.close()

            if options.bom_only:
                exit(0)

            for i, bom_row in enumerate(bom_table_bot):
                print("Plotting bottom page (%d/%d) %s / %s" % (i+1, len(bom_table_bot),
                                                         bom_row[2], bom_row[1]))
                create_board_figure(pcb, bom_row, options.boards, layer=pcbnew.B_Cu)
                pdf.savefig()
                plt.close()

            for i, bom_row in enumerate(bom_table_top):
                print("Plotting top page (%d/%d) %s / %s" % (i+1, len(bom_table_top),
                                                         bom_row[2], bom_row[1]))
                create_board_figure(pcb, bom_row, options.boards, layer=pcbnew.F_Cu)
                pdf.savefig()
                plt.close()
            print("Output written to %s" % fname_out)
    else:
        if options.fsort:
            shutil.rmtree('smt', ignore_errors=True)
            shutil.rmtree('tht', ignore_errors=True)
            os.makedirs('smt/small')
            os.makedirs('tht/connectors')

        for i, bom_row in enumerate(bom_table_top):
            if bom_row[3] == 1:
                mtype = 'smt'
            else:
                mtype = 'tht'

            msubtype = ''
            for r in bom_row[4]:
                if re.match('^[X][0-9]', r) and mtype == 'tht':
                    msubtype = 'connectors/'
            if re.match('^[RCD][0-9]{4}$', bom_row[2]) or re.match('^D_SOD[0-9]', bom_row[2]) or re.match('^[RC]NET', bom_row[2]):
                msubtype = 'small/'

            if options.fsort:
                fname_out = mtype + "/" + msubtype + bom_row[2].replace('/', '_') + '_' + bom_row[1].replace(
                    '/', '_') + "_" + os.path.splitext(file)[0] + "_picknplace.pdf"
            else:
                fname_out = bom_row[2].replace('/', '_') + '_' + bom_row[1].replace(
                    '/', '_') + "_" + mtype + "_" + os.path.splitext(file)[0] + "_picknplace.pdf"

            print("Plotting (%d/%d) %s / %s to %s" %
                  (i+1, len(bom_table_top), bom_row[2], bom_row[1], fname_out))
            with PdfPages(fname_out) as pdf:
                create_board_figure(pcb, bom_row, layer=pcbnew.F_Cu)
                pdf.savefig()
                plt.close()
