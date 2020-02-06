import math
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
    


def readColorFile(filename):
    black = "0.0  0.0  0.0"
    setcolor = {'A':black, 'C':black, 'G':black, 'T':black, 'a':black, 'c':black, 'g':black, 't':black}
    f = open(colorfile, 'r')    
    lines = f.readlines()    
    for line in lines:
        parts = line.strip().split(' ')
        setcolor[parts[0]] = parts[1] + "  " +  parts[2] + "  " +  parts[3] + "  setrgbcolor\n"
    f.close()
    return setcolor
    
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


def pt2mm(pt):
    return pt * 0.352778

def s2ps2(
        sequenceData,                  # dictionary {"sequence":<string>, "name":<name of organism>}
        filename,                       # name of the file to write (e.g., "myseq.ps")
        colorfile       = "colors.txt",
          
        paper           = "Letter",     # Type of paper (used to determine paper size)
        d_mm            = 1.0,          # Size in mm of the squares

        nameFontsize    = 16,

        nonCoding       = False,        # If True then noncoding regions are colored differently
        legend_width_cm = 3.2,          # Width of the bottom color legend
        frameboxes      = 2,            # number of frames to put around the graphic
        gap_mm          = 1.0,          # size of gap to put between frames
        framestroke     = 0.5           # stroke width of the frames
        ):

        seq      = sequence['sequence']
        name     = sequence['name']
        setcolor = readColorFile(colorfile)

        f = open(filename + ".ps", 'w')
        paper_width_pt, paper_height_pt = paper2pt(paper)
        paper_width_mm  = pt2mm(paper_width_pt)
        paper_height_mm = pt2mm(paper_height_pt)

        width_mm = paper_width_mm - 2 * left_margin_mm
        height_mm = paper_height_mm - 2*top_margin_mm

        cols = math.floor(width_mm / d_mm)
        rows = math.ceil(len(sequence) / cols)
        last_row      = len(sequence) - rows * cols

        rows_per_page = math.floor(height_mm / d_mm)
        pages         = len(seq) / rows_per_page
        
    f.write("%!\ n")
    f.write("% DNA sequence visualization\n\n")
    f.write("/mm {2.83464  mul} def\n\n")  # Convert mm to pt

    left_margin_mm = 0.5 * paper_width_mm - 0.5 * width_mm
    bottom_margin_mm = (paper_height_mm - rows*d_mm) / 2
    
    bp = 0         
         pass