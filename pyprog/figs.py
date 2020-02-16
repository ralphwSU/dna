import math
import re
import genelink as g

# width  = width of image in cm
# height = height of image in cm
# d      = dimension of square
def s2t(sequence, filename, width_cm, d_mm):
    f = open(filename + ".tex", 'w')
    width_mm = width_cm * 10
    cols = math.floor(width_mm / d_mm)
    rows = math.floor(len(sequence) / cols)
    last_row      = len(sequence) - rows * cols

    f.write("\\documentclass[12pt]{article}\n")
    f.write("\\usepackage[utf8]{inputenc}\n")
    f.write("\setlength{\\textwidth}{6.4in}\n")
    f.write("\\setlength{\\textheight}{9.0in}\n")
    f.write("\\setlength{\\oddsidemargin}{0pt}\n")
    f.write("\\setlength{\\evensidemargin}{0pt}\n")
    f.write("\\setlength{\\topmargin}{0pt}\n")
    f.write("\\setlength{\\hoffset}{.05in}\n")
    f.write("\\setlength{\\voffset}{-.5in}\n")
    f.write("\\setlength{\\parskip}{5pt}\n")
    f.write("\\setlength{\\parindent}{0pt}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{color}\n")
    f.write("\\usepackage{graphicx}\n")     
    f.write("\\usepackage[T1]{fontenc}\n")
    f.write("\\usepackage{ifthen}\n")
    f.write("\\usepackage{pstricks,pst-grad,pst-text,pst-node,pst-func,pst-plot}\n")
    f.write("\\definecolor{A}{rgb}{1,0,0}\n")
    f.write("\\definecolor{C}{rgb}{0,1,0}\n")
    f.write("\\definecolor{G}{rgb}{0,0,1}\n")
    f.write("\\definecolor{T}{rgb}{.5,.5,.5}\n")        
    f.write("\\begin{document}")
    f.write("\\pagestyle{empty}")

    print(str(rows) +  ", " +  str(cols))
    
    bp = 0
    f.write("\\begin{center}\n")
    f.write("\\begin{pspicture}(0,0)(" + str(cols*d_mm / 10) + "," + str(rows*d_mm / 10) + ")\n")
    for row in range(rows):
        for col in range(cols):
            s = "\psframe*[linecolor=" + sequence[bp] + "]("
            s = s + str((col-1)*d_mm / 10)
            s = s + ","
            s = s + str((row-1)*d_mm / 10)
            s = s + ")("
            s = s + str(col*d_mm / 10)
            s = s + ","
            s = s + str(row*d_mm / 10)
            s = s + ")"
            f.write(s)
            f.write("\n")
            bp = bp + 1
    f.write("\\end{pspicture}")    
    f.write("\\end{center}")    
    f.write("\\vfill\\eject")
    f.write("\\end{document}")

def F(x):
    return '{0:0.5f}'.format(x)

def F3(x):
    return '{0:0.3f}'.format(x)

# Convert gene sequence to a postscript file
def s2ps(sequence, sequenceName, filename, width_cm, d_mm, colorfile = "colors.txt", fontsize = 16):
    # Read colors
    black = "0.0  0.0  0.0"
    setcolor = {'A':black, 'C':black, 'G':black, 'T':black, 'a':black, 'c':black, 'g':black, 't':black}
    
    f = open(colorfile, 'r')    
    lines = f.readlines()    
    for line in lines:
        parts = line.strip().split(' ')
        setcolor[parts[0]] = parts[1] + "  " +  parts[2] + "  " +  parts[3] + "  setrgbcolor\n"
    f.close()

    f = open(filename + ".ps", 'w')
    paper_width_mm = 8.5 * 2.54 * 10
    paper_height_mm = 11 * 2.54 * 10
    width_mm = width_cm * 10
    cols = math.floor(width_mm / d_mm)
    rows = math.ceil(len(sequence) / cols)
    last_row      = len(sequence) - rows * cols

    f.write("%!\ n")
    f.write("% DNA sequence visualization\n\n")
    f.write("/mm {2.83464  mul} def\n\n")  # Convert mm to pt

    left_margin_mm = 0.5 * paper_width_mm - 0.5 * width_mm
    bottom_margin_mm = (paper_height_mm - rows*d_mm) / 2
    
    bp = 0
    for row in range(rows - 1): # Do the last row separately since it may not be filled
        for col in range(cols):
           f.write(setcolor[sequence[bp]])
           f.write("newpath\n")
           left   = left_margin_mm + col * d_mm
           right  = left_margin_mm + (col + 1) * d_mm
           top    = bottom_margin_mm + (rows - row) * d_mm
           bottom = bottom_margin_mm + (rows - row - 1) * d_mm
           f.write(F(left)  + " mm " + '{0:.5f}'.format(top)    + " mm moveto\n")
           f.write('{0:.5f}'.format(right) + " mm " + '{0:.5f}'.format(top)    + " mm lineto\n")
           f.write('{0:.5f}'.format(right) + " mm " + '{0:.5f}'.format(bottom) + " mm lineto\n")
           f.write('{0:.5f}'.format(left)  + " mm " + '{0:.5f}'.format(bottom) + " mm lineto\n")           
           f.write("closepath\n")
           f.write("fill\n")
           bp = bp + 1
    remaining_pairs = len(sequence) - bp
    row = rows - 1 
    for col in range(remaining_pairs):
           f.write(setcolor[sequence[bp]])
           f.write("newpath\n")
           left   = left_margin_mm + col * d_mm
           right  = left_margin_mm + (col + 1) * d_mm
           top    = bottom_margin_mm + (rows - row) * d_mm
           bottom = bottom_margin_mm + (rows - row - 1) * d_mm
           f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
           f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
           f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
           f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")           
           f.write("closepath\n")
           f.write("fill\n")       
           bp = bp + 1

    # Print two boxes around the image
    gap_mm = 1.0 
    f.write("% Print a box\n")
    f.write("0.0  0.0  0.0  setrgbcolor\n")
    f.write("0.5 setlinewidth\n")
    f.write("newpath\n")
    f.write(F(left_margin_mm)                  + " mm " + F(bottom_margin_mm) + " mm moveto\n")
    f.write(F(left_margin_mm)                  + " mm " + F(paper_height_mm - bottom_margin_mm) + " mm lineto\n")
    f.write(F(paper_width_mm - left_margin_mm) + " mm " + F(paper_height_mm - bottom_margin_mm) + " mm lineto\n")
    f.write(F(paper_width_mm - left_margin_mm) + " mm " + F(bottom_margin_mm) + " mm lineto\n")
    f.write("closepath\n")
    f.write("stroke\n")
    #
    f.write("newpath\n")
    f.write(F(left_margin_mm - gap_mm) + " mm " + F(bottom_margin_mm - gap_mm) + " mm moveto\n")
    f.write(F(left_margin_mm - gap_mm) + " mm " + F(paper_height_mm - bottom_margin_mm + gap_mm) + " mm lineto\n")
    f.write(F(paper_width_mm - left_margin_mm + gap_mm) + " mm " + F(paper_height_mm - bottom_margin_mm + gap_mm) + " mm lineto\n")
    f.write(F(paper_width_mm - left_margin_mm + gap_mm) + " mm " + F(bottom_margin_mm - gap_mm) + " mm lineto\n")
    f.write("closepath\n")
    f.write("stroke\n")

    # Print the name of the organism
    f.write("0.0  0.0  0.0  setrgbcolor\n")
    f.write("/Times-Roman findfont\n")
    f.write(str(fontsize) + " scalefont\n")
    f.write("setfont\n")
    bottom_of_text = bottom_margin_mm - 2.5 * fontsize * (2.54 * 10 / (2 * 72)) 
    f.write(F(left_margin_mm) + " mm " + F(bottom_of_text) + " mm moveto\n")
    f.write("(" + sequenceName + ") show\n")

    f.write("0.5  0.5  0.5  setrgbcolor\n")    
    h = 1.5 * fontsize * 2.54 * 10 / 72  # distance of color box below base pair letter    
    f.write(F(left_margin_mm) + " mm " + F(bottom_of_text - h) + " mm moveto\n")
    f.write("(" + str(len(sequence)) + " base pairs) show\n")    
    
    # Now print the legend
    legend_right_offset_mm = fontsize * 2.54 * 10 / 72
    w = 0.68 * fontsize * 2.54 * 10 / 72  # width of legend color box
    h = 0.28 * fontsize * 2.54 * 10 / 72  # distance of color box below base pair letter
    legend_width_mm        = 4.0 * 8
    f.write("0.3  0.3  0.3  setrgbcolor\n")
    f.write("/Times-Roman findfont\n")
    f.write(str(fontsize) + " scalefont\n")
    f.write("setfont\n")
    bottom_of_text = bottom_margin_mm - 2.5 * fontsize * (2.54 * 10 / (2 * 72))

    x_pos_mm = paper_width_mm - left_margin_mm - legend_right_offset_mm - 3 * legend_width_mm / 4
    f.write(F(x_pos_mm) + " mm " + F(bottom_of_text) + " mm moveto\n")
    f.write("(A) show\n")
    f.write(setcolor['A'])
    f.write("newpath\n")
    left   = x_pos_mm
    right  = x_pos_mm + w
    top    = bottom_of_text - h
    bottom = bottom_of_text - h - w
    f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
    f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
    f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
    f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")           
    f.write("closepath\n")
    f.write("fill\n")

    f.write("0.3  0.3  0.3  setrgbcolor\n")    
    x_pos_mm = paper_width_mm - left_margin_mm - legend_right_offset_mm - 2 * legend_width_mm / 4
    f.write(F(x_pos_mm) + " mm " + F(bottom_of_text) + " mm moveto\n")
    f.write("(C) show\n")
    f.write(setcolor['C'])    
    f.write("newpath\n")
    left   = x_pos_mm
    right  = x_pos_mm + w
    top    = bottom_of_text - h
    bottom = bottom_of_text - h - w
    f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
    f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
    f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
    f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")           
    f.write("closepath\n")
    f.write("fill\n")

    f.write("0.3  0.3  0.3  setrgbcolor\n")        
    x_pos_mm = paper_width_mm - left_margin_mm - legend_right_offset_mm -  legend_width_mm / 4    
    f.write(F(x_pos_mm) + " mm " + F(bottom_of_text) + " mm moveto\n")
    f.write("(G) show\n")
    f.write(setcolor['G'])    
    f.write("newpath\n")
    left   = x_pos_mm
    right  = x_pos_mm + w
    top    = bottom_of_text - h
    bottom = bottom_of_text - h - w
    f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
    f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
    f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
    f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")           
    f.write("closepath\n")
    f.write("fill\n")

    f.write("0.3  0.3  0.3  setrgbcolor\n")        
    x_pos_mm = paper_width_mm - left_margin_mm - legend_right_offset_mm
    f.write(F(x_pos_mm) + " mm " + F(bottom_of_text) + " mm moveto\n")
    f.write("(T) show\n")            
    f.write(setcolor['T'])    
    f.write("newpath\n")
    left   = x_pos_mm
    right  = x_pos_mm + w
    top    = bottom_of_text - h
    bottom = bottom_of_text - h - w
    f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
    f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
    f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
    f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")           
    f.write("closepath\n")
    f.write("fill\n")
    
    f.write("showpage")
    f.close()
    


def readColorFile(colorfile):
    black = "0.0  0.0  0.0"
    setcolor = {'A':black, 'C':black, 'G':black, 'T':black, 'a':black, 'c':black, 'g':black, 't':black}
    f = open(colorfile, 'r')    
    lines = f.readlines()    
    for line in lines:
        parts = line.strip().split(' ')
        setcolor[parts[0]] = parts[1] + "  " +  parts[2] + "  " +  parts[3] + "  setrgbcolor\n"
    f.close()
    return setcolor

# Compute the paper size in points from the name of the paper size
def paper2pt(paper = "Letter"):
    if paper == "Comm_10_Envelope":
       return [297, 684]
    if paper == "C5_Envelope":
        return [461, 468]
    if paper == "DL_Envelope":
        return [312, 624]
    if paper == "Folio":
        return [595, 935]
    if paper == "Executive":
        return [522, 756]
    if paper == "Letter":
        return [612, 792]
    if paper == "Legal":
        return [612, 1008]
    if paper == "Ledger":
        return [1224, 792]
    if paper == "Tabloid":
        return [792, 1224]
    if paper == "A0":
        return [2384, 3370]
    if paper == "A1":
        return [1684, 2384]
    if paper == "A2":
        return [1191, 1684]
    if paper == "A3":
        return [842, 1191]
    if paper == "A4":
        return [595, 842]
    if paper == "A5":
        return [420, 595]
    if paper == "A6":
        return [297, 420]
    if paper == "A7":
        return [210, 297]
    if paper == "A8":
        return [148, 210]
    if paper == "A9":
        return [105, 148]
    if paper == "B0":
        return [105, 148]
    if paper == "B1":
        return [2064, 2920]
    if paper == "B2":
        return [1460, 2064]
    if paper == "B3":
        return [1032, 1460]
    if paper == "B4":
        return [729, 1032]
    if paper == "B5":
        return [516, 729]
    if paper == "B6":
        return [363, 516]
    if paper == "B7":
        return [258, 363]
    if paper == "B8":
        return [181, 258]
    if paper == "B9":
        return [127, 181]
    if paper == "B10":
        return [91, 127]


# Convert points to mm.    
def pt2mm(pt):
    return pt * 0.352777778

# Given the name of a sequence file, generate the ps image file.
def ps(
        fastafile,
        colorfile         = "colorsBlue3.txt",
          
        paper             = "Letter",     # Type of paper (used to determine paper size)
        figure_width_cm   = 13.0,         # Width of the genome plot.
        top_margin_cm     = 1.5,          # Width of the top margin on each page.
        legend_width_cm   = 3.2,          # Width of the bottom color legend        
        d_mm              = 1.0,          # Size in mm of the squares

        fontName          = "Times-Roman",
        nameFontSize      = 16,
        nameShiftFactor   = 1.5,          # Controls how much to shift the organism name above the image
        pageFontSize      = 10,           # Size of the page number font
        pageShiftFactor_h = 0.25,         # Controls how far into the right margin to print page number
        pageShiftFactor_v = 0.6,         # Controls how far below the baseline to print the page number
        statFontSize      = 12,           # Size of the font used to print statistics
        statShiftFactor_v = 0.3,          # Controls how far away from the frame to shift the stats
        statShiftFactor_h = 0.1,          # Controls how far away from the frame to shift the stats 
        frameSepFactor    = 0.5,
        
        nonCoding         = False,        # If True then noncoding regions are colored differently

        frameboxes        = 2,            # number of frames to put around the graphic
        gap_mm            = 1.0,          # size of gap to put between frames
        framestroke       = 0.5           # stroke width of the frames
        ):

   s = g.readseq2(fastafile)

   filename = re.sub(' +', '_', s['name'])
   
   s2ps2(s,
         filename,
         colorfile,
         paper,
         figure_width_cm,
         top_margin_cm,
         legend_width_cm,
         d_mm,
         fontName,
         nameFontSize,
         nameShiftFactor,
         pageFontSize,
         pageShiftFactor_h,
         pageShiftFactor_v,
         statFontSize,
         statShiftFactor_v,
         statShiftFactor_h,
         frameSepFactor,
         nonCoding,
         frameboxes,
         gap_mm,
         framestroke)


def write_centering_macro(f):
   f.write("% Macro for centering text.\n")
   f.write("/ceshow { % (string) fontsize fontname x y\n")
   f.write("gsave\n")
   f.write("moveto % s fs fn\n")
   f.write("findfont exch scalefont setfont % s\n")
   f.write("gsave\n")
   f.write("dup false charpath flattenpath\n")
   f.write("pathbbox % s x0 y0 x1 y1\n")
   f.write("grestore\n")
   f.write("3 -1 roll sub % s x0 x1 dy\n")
   f.write("3 1 roll sub % s dy -dx\n")
   f.write("2 div exch % s -dx/2 dy\n")
   f.write("-2 div % s -dx/2 -dy\n")
   f.write("rmoveto show\n")
   f.write("grestore\n")
   f.write("} bind def\n\n")

def write_center_align_macro(f):
   f.write("/cshow { dup stringwidth pop 2 div neg 0 rmoveto show } def\n\n")
   
def write_right_align_macro(f):
   f.write("/rshow { dup stringwidth pop neg 0 rmoveto show } def\n\n")
   
def write_rotate_text_macro(f, fontName):
   f.write("% Macro for typesetting rotated text\n")
   f.write("/outputtext {\n")
   f.write("/data exch def\n")
   f.write("/rot exch def\n")
   f.write("/xfont exch def\n")
   f.write("/y1 exch def\n")
   f.write("/x1 exch def\n")
   f.write("/" + fontName + " findfont\n")
   f.write("xfont scalefont\n")
   f.write("setfont\n")
   f.write("x1 y1 moveto\n")
   f.write("rot rotate\n")
   f.write("data show\n")
   f.write("rot neg rotate\n")
   f.write("} def\n\n")
   
# Convert a genome to a postscript file
def s2ps2(
        sequenceData,                   # dictionary {"sequence":<string>, "name":<name of organism>}
        filename,                       # name of the file to write (e.g., "myseq.ps")
        colorfile         = "colorsBlue3.txt",
          
        paper             = "Letter",     # Type of paper (used to determine paper size)
        figure_width_cm   = 13.0,         # Width of the genome plot.
        top_margin_cm     = 1.5,          # Width of the top margin on each page.
        legend_width_cm   = 3.2,          # Width of the bottom color legend        
        d_mm              = 1.0,          # Size in mm of the squares

        fontName          = "Times-Roman",
        nameFontSize      = 16,           # Size of the font used for the organism name
        nameShiftFactor   = 1.5,          # Controls how much to shift the organism name above the image
        pageFontSize      = 10,           # Size of the page number font
        pageShiftFactor_h = 0.25,         # Controls how far into the right margin to print page number
        pageShiftFactor_v = 0.6,          # Controls how far below the baseline to print the page number

        statFontSize      = 12,           # Size of the font used to print statistics
        statShiftFactor_v = 0.3,          # Controls how far away from the frame to shift the stats
        statShiftFactor_h = 0.1,          # Controls how far away from the frame to shift the stats
        frameSepFactor    = 0.5,          # Controls separation between the two frames

        nonCoding         = False,        # If True then noncoding regions are colored differently

        frameboxes        = 2,            # number of frames to put around the graphic
        gap_mm            = 1.0,          # size of gap to put between frames
        framestroke       = 0.5           # stroke width of the frames
        ):

    # First set a default color then read the color file.
    black = "0.0  0.0  0.0"
    setcolor = {'A':black, 'C':black, 'G':black, 'T':black, 'a':black, 'c':black, 'g':black, 't':black}
    f = open(colorfile, 'r')    
    lines = f.readlines()    
    for line in lines:
        parts = line.strip().split(' ')
        setcolor[parts[0]] = parts[1] + "  " +  parts[2] + "  " +  parts[3] + "  setrgbcolor\n"
    f.close()

    sequence  = sequenceData['sequence']
    sequence1 = g.seq2c2(sequenceData)
    stats     = g.stat(sequence1)
    sequence = sequence1['cs']
    name     = sequenceData['name']
    setcolor = readColorFile(colorfile)

    paper_width_pt, paper_height_pt = paper2pt(paper)
    paper_width_mm   = pt2mm(paper_width_pt)
    paper_height_mm  = pt2mm(paper_height_pt)

    # What if figure_width_mm > paper_width_mm?
    # How to determine margin sizes?
    figure_width_mm  = figure_width_cm * 10    
    horiz_margin_mm  = (paper_width_mm - figure_width_mm) / 2.0
    figure_left_mm   = horiz_margin_mm
    figure_right_mm  = figure_left_mm + figure_width_mm
    
    top_margin_mm    = top_margin_cm * 10
    figure_top_mm    = paper_height_mm - top_margin_mm
    

    cols = math.floor(figure_width_mm / d_mm)
    figure_right_mm  = figure_left_mm + cols * d_mm
    figure_width_mm  = figure_right_mm - figure_left_mm
    
    rows = math.ceil(len(sequence) / cols)
    bp_on_last_row      = len(sequence) - (rows - 1) * cols
    empty_bp_on_last_row = cols - bp_on_last_row

    figure_height_mm  = paper_height_mm - 2 * top_margin_mm
    rows_per_page = math.floor(figure_height_mm / d_mm)
    
    figure_bot_mm     = figure_top_mm - rows_per_page * d_mm
    figure_height_mm  = figure_top_mm - figure_bot_mm    

    bp_per_page   = cols * rows_per_page
    pages         = math.ceil(len(sequence) / bp_per_page)
    bp_on_last_page = len(sequence) - bp_per_page * (pages - 1)
    rows_on_last_page = math.ceil(bp_on_last_page / cols)
    
    figure_top_mm    = paper_height_mm - top_margin_mm
    name_y_mm        = figure_top_mm + nameShiftFactor * pt2mm(nameFontSize)

    print("paper_width_mm       = " + str(paper_width_mm))
    print("paper_height_mm      = " + str(paper_height_mm))
    print("cols                 = " + str(cols))    
    print("rows                 = " + str(rows))
    print("rows_per_page        = " + str(rows_per_page))
    print("bp_per_page          = " + str(bp_per_page))    
    print("pages                = " + str(pages))
    print("bp_on_last_page      = " + str(bp_on_last_page))
    print("bp_on_last_row       = " + str(bp_on_last_row))
    print("rows_on_last_page    = " + str(rows_on_last_page))
    print("empty_bp_on_last_row = " + str(empty_bp_on_last_row))
    print("figure_width_mm      = " + str(figure_width_mm))
    print("figure_left_mm       = " + str(figure_left_mm))
    print("figure_right_mm      = " + str(figure_right_mm))
    print("figure_height_mm     = " + str(figure_height_mm))    
    print("figure_top_mm        = " + str(figure_top_mm))
    print("figure_bot_mm        = " + str(figure_bot_mm))        
    
    f = open(filename + ".ps", 'w')    
    f.write("%!\n")
    f.write("% DNA sequence visualization for " + name + "\n\n")

    # Write postscript macros
    f.write("% Macro to convert mm to pt.\n")
    f.write("/mm {2.83464  mul} def\n\n")  
    write_centering_macro(f)
    write_rotate_text_macro(f, fontName)
    write_right_align_macro(f)
    write_center_align_macro(f)    
    

    # Make the legend    

    bp = 0
    for page in range(1, pages + 1):
       f.write("%%Page: " + str(page) + " " + str(pages) + "\n\n")

       # Write the organism name       
       f.write("% Write the organism name.\n")
       # Set color of the organism name to black on the first
       # page and gray on subsequent ones.
       if page == 1:
          f.write("0.0  0.0  0.0  setrgbcolor\n")
       else:
          f.write("0.6  0.6  0.6  setrgbcolor\n")
       f.write("/" + fontName + " findfont\n")
       f.write(str(nameFontSize) + " scalefont\n")
       f.write("setfont\n")

       f.write("(" + name + ") ")
       f.write(str(nameFontSize) + " /" + fontName + " ")
       f.write(str(paper_width_mm / 2.0) + " mm ")
       f.write(str(name_y_mm) + " mm ")
       f.write("ceshow\n\n")


       rows_this_page = rows_per_page
       if page == pages:
           rows_this_page = rows_on_last_page
         
       for row in range(rows_this_page):
          cols_this_row = cols
          if (page == pages) and (row == (rows_on_last_page - 1)):
              cols_this_row = bp_on_last_row
          for col in range(cols_this_row):
             f.write(setcolor[sequence[bp]])
             f.write("newpath\n")
             left   = figure_left_mm + col * d_mm
             right  = figure_left_mm + (col + 1) * d_mm

             top    = figure_top_mm  - row * d_mm
             bottom = figure_top_mm  - (row + 1) * d_mm             
             
             #top    = figure_bot_mm  + (rows_this_page - row) * d_mm
             #bottom = figure_bot_mm  + (rows_this_page - row -1) * d_mm
             f.write(F(left)  + " mm " + F(top)    + " mm moveto\n")
             f.write(F(right) + " mm " + F(top)    + " mm lineto\n")
             f.write(F(right) + " mm " + F(bottom) + " mm lineto\n")
             f.write(F(left)  + " mm " + F(bottom) + " mm lineto\n")  
             f.write("closepath\n")
             f.write("fill\n")
             bp = bp + 1
       f.write("\n")
             
       # Make the frame(s)
       f.write("% Print frame(s) around the image.\n")
       f.write("0.0  0.0  0.0  setrgbcolor\n")
       f.write("0.5 setlinewidth\n")
       f.write("newpath\n")
       # bottom should change on last page
       if page == pages:
           figure_bot_mm = figure_top_mm - rows_on_last_page * d_mm
       f.write(F(figure_left_mm)  + " mm " + F(figure_bot_mm) + " mm moveto\n")
       f.write(F(figure_right_mm) + " mm " + F(figure_bot_mm) + " mm lineto\n")
       f.write(F(figure_right_mm) + " mm " + F(figure_top_mm) + " mm lineto\n")
       f.write(F(figure_left_mm)  + " mm " + F(figure_top_mm) + " mm lineto\n")
       f.write("closepath\n")
       f.write("stroke\n\n")       

       f.write("0.0  0.0  0.0  setrgbcolor\n")
       f.write("0.5 setlinewidth\n")
       f.write("newpath\n")
       # bottom should change on last page
       if page == pages:
           figure_bot_mm = figure_top_mm - rows_on_last_page * d_mm
       x_left  = figure_left_mm  - frameSepFactor * d_mm
       x_right = figure_right_mm + frameSepFactor * d_mm
       y_top   = figure_top_mm   + frameSepFactor * d_mm
       y_bot   = figure_bot_mm   - frameSepFactor * d_mm       
       f.write(F(x_left)  + " mm " + F(y_bot) + " mm moveto\n")
       f.write(F(x_right) + " mm " + F(y_bot) + " mm lineto\n")
       f.write(F(x_right) + " mm " + F(y_top) + " mm lineto\n")
       f.write(F(x_left)  + " mm " + F(y_top) + " mm lineto\n")
       f.write("closepath\n")
       f.write("stroke\n\n")       
       
       f.write("% Print number of base pairs.\n")
       f.write("0.5  0.5  0.5  setrgbcolor\n")
       f.write("/" + str(fontName) + " findfont\n")
       f.write(str(statFontSize) + " scalefont\n")
       f.write("setfont\n")       
       x_pos = horiz_margin_mm
       y_pos = top_margin_mm - statShiftFactor_v * top_margin_mm # should depend on font size
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm moveto\n")
       f.write("(" + str(len(sequence)) + " base pairs) show\n\n ")    
       
       # Print  "print page / pages" at the lower left corner      
       f.write("% Print page / pages at the lower left corner\n")
       f.write("0.4  0.4  0.4  setrgbcolor\n")
       f.write("/" + str(fontName) + " findfont\n")
       f.write(str(pageFontSize) + " scalefont\n")
       f.write("setfont\n")
       text = "page " + str(page) + "/" + str(pages)
       x_pos = (paper_width_mm - horiz_margin_mm) + pageShiftFactor_h * horiz_margin_mm
       y_pos = top_margin_mm  - pageShiftFactor_v * top_margin_mm       
       f.write(str(x_pos) + " mm " + str(y_pos) + " mm moveto\n")
       f.write("(" + text + ") show\n\n")

       # Write statistics
       x_pos = figure_left_mm - statShiftFactor_h * horiz_margin_mm
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm
       f.write("% Print percent of each coding/non-coding base molecule\n")
       f.write("0.3  0.3  0.3  setrgbcolor\n")
       text = "(A: " + F3(stats['A']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + figure_height_mm / 8
       text = "(C: " + F3(stats['C']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 2 * figure_height_mm / 8
       text = "(G: " + F3(stats['G']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 3 * figure_height_mm / 8
       text = "(T: " + F3(stats['T']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 4 * figure_height_mm / 8
       text = "(a: " + F3(stats['a']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 5 * figure_height_mm / 8
       text = "(c: " + F3(stats['c']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 6 * figure_height_mm / 8
       text = "(g: " + F3(stats['g']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")

       y_pos = figure_bot_mm + 7 * figure_height_mm / 8
       text = "(t: " + F3(stats['t']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n\n")       

       # Write the color legend
       f.write("% Write the color legend\n")
       legendBoxSize_mm = 3
       f.write(setcolor['A'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n") 

       f.write(setcolor['C'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")

       f.write(setcolor['G'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + 2 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")        

       f.write(setcolor['T'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + 3 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")

       f.write(setcolor['a'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm + 4 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n") 

       f.write(setcolor['c'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + 5 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")

       f.write(setcolor['g'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + 6 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")        

       f.write(setcolor['t'])
       f.write("newpath\n")
       x_pos = figure_left_mm - 2.5 * statShiftFactor_h * horiz_margin_mm       
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm  + 7 * figure_height_mm / 8
       f.write(F(x_pos)                    + " mm " + F(y_pos) + " mm moveto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos) + " mm lineto\n")
       f.write(F(x_pos - legendBoxSize_mm) + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")
       f.write(F(x_pos)                    + " mm " + F(y_pos + legendBoxSize_mm) + " mm lineto\n")       
       f.write("closepath\n")
       f.write("fill\n")        
       
       # Write the coding fraction
       f.write("% Write coding and noncoding fractions\n")
       f.write("0.3  0.3  0.3  setrgbcolor\n")
       f.write("/" + str(fontName) + " findfont\n")
       f.write(str(statFontSize) + " scalefont\n")
       f.write("setfont\n")       
       x_pos = paper_width_mm / 2
       y_pos = top_margin_mm - statShiftFactor_v * top_margin_mm 
       text = "(coding : " + F3(stats['i']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm moveto\n")
       f.write(text + " cshow\n")
       
       # Right justify the non-coding string
       x_pos = figure_right_mm
       text = "(noncoding : " + F3(stats['e']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm moveto\n")
       f.write(text + " rshow\n\n")

       # Write the amino acid percentages
       x_pos = figure_right_mm + 1.5 * statShiftFactor_h * horiz_margin_mm
       figure_bot_mm = figure_top_mm - rows_per_page * d_mm 
       y_pos = figure_bot_mm
       f.write("% Print percent of each amino acid in the coding section\n")
       # 1
       f.write("0.3  0.3  0.3  setrgbcolor\n")
       text = "(A: " + F3(stats['A1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 2
       y_pos = figure_bot_mm + figure_height_mm / 10
       text = "(C: " + F3(stats['C1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 3       
       y_pos = figure_bot_mm + 2 * figure_height_mm / 10
       text = "(D: " + F3(stats['D1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 4       
       y_pos = figure_bot_mm + 3 * figure_height_mm / 10
       text = "(E: " + F3(stats['E1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 5       
       y_pos = figure_bot_mm + 4 * figure_height_mm / 10
       text = "(F: " + F3(stats['F1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 6       
       y_pos = figure_bot_mm + 5 * figure_height_mm / 10
       text = "(G: " + F3(stats['G1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 7       
       y_pos = figure_bot_mm + 6 * figure_height_mm / 10
       text = "(H: " + F3(stats['H1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 8
       y_pos = figure_bot_mm + 7 * figure_height_mm / 10
       text = "(I: " + F3(stats['I1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 9       
       y_pos = figure_bot_mm + 8 * figure_height_mm / 10
       text = "(K: " + F3(stats['K1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 10       
       y_pos = figure_bot_mm + 9 * figure_height_mm / 10
       text = "(L: " + F3(stats['L1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")       

       # 11
       x_pos = x_pos + 2 * statShiftFactor_h * horiz_margin_mm
       y_pos = figure_bot_mm #+ figure_height_mm 
       text = "(M: " + F3(stats['M1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 12      
       y_pos = figure_bot_mm + 1 * figure_height_mm / 10
       text = "(N: " + F3(stats['N1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 13       
       y_pos = figure_bot_mm + 2 * figure_height_mm / 10
       text = "(P: " + F3(stats['P1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 14       
       y_pos = figure_bot_mm + 3 * figure_height_mm / 10
       text = "(Q: " + F3(stats['Q1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 15       
       y_pos = figure_bot_mm + 4 * figure_height_mm / 10
       text = "(R: " + F3(stats['R1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 16       
       y_pos = figure_bot_mm + 5 * figure_height_mm / 10
       text = "(S: " + F3(stats['S1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 17
       y_pos = figure_bot_mm + 6 * figure_height_mm / 10
       text = "(T: " + F3(stats['T1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 18      
       y_pos = figure_bot_mm + 7 * figure_height_mm / 10
       text = "(V: " + F3(stats['V1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")
       # 19       
       y_pos = figure_bot_mm + 8 * figure_height_mm / 10
       text = "(Y: " + F3(stats['Y1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")       
       # 20
       y_pos = figure_bot_mm + 9 * figure_height_mm / 10
       text = "(W: " + F3(stats['W1']) + ")"
       f.write(F(x_pos) + " mm " + F(y_pos) + " mm ")
       f.write(str(statFontSize) + " 90 " + text + " outputtext\n")       
       
       f.write("showpage\n\n")
    f.close()

import os
dataDir = "../data/"
figDir  = "../figs/"
import glob
def makebook():
    fasta_files = []
    for file in glob.glob(dataDir + "*.fasta"):
       fasta_files.append(file)
    for file in fasta_files:
        fig_file = figDir + file[8:len(file)]
        fig_file = fig_file[0:(len(fig_file) - 6)]
        fig_file = fig_file
        s = g.readseq2(file)        
        s2ps2(s, filename = fig_file, d_mm = 1.2, colorfile='colorsBlue2.txt', top_margin_cm = 2.0)
        os.system("ps2pdf " + fig_file + ".ps")
        fig_file = fig_file[8:len(fig_file)]
        os.system("mv " + fig_file + ".pdf " + figDir)
    os.chdir(figDir)
    pdf_files = []
    for file in glob.glob("*.pdf"): #figDir + "*.pdf"):
       pdf_files.append(file)
    print(pdf_files)
    
    os.system("pdftk A=" + pdf_files[0] + " B=" + pdf_files[1] + " cat A B output tmp.pdf")
    for i in range(2, len(pdf_files)):
        os.system("pdftk A=tmp.pdf B=" + pdf_files[i] + " cat A B output tmp.pdf")
    os.chdir("../pyprog")        


