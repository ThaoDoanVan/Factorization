#
# apt-get install python-tk
# python3 projectUI.py
#
import os
from threading import Thread
from ctypes import *
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
from tkinter.filedialog import askopenfile

def select_file():
	global linkFile
	global inputText
	filename = askopenfile(mode='r', filetypes=[('Data files', '*.txt')])
	if filename is not None:
		linkFile = filename.name
		inputFilenameText.set("Filename: "+os.path.split(linkFile)[1])
		linkFile = os.path.split(linkFile)[1]
		startButton.config(state='active')
		with open(linkFile, 'r') as file:
			data = file.read().splitlines()
		for i in range(5):
			if i < len(data):
				inputText[i].set(data[i])
			else:
				inputText[i].set("")

def changeSampleInput():
	if chooseMethod.get() == 0:
		sampleInputText1.set("n = 2047")
		sampleInputText2.set("p_max = 113")
		sampleInputText3.set("")
		sampleInputText4.set("")
		sampleInputText5.set("")
	if chooseMethod.get() == 1:
		sampleInputText1.set("n = 53110")
		sampleInputText2.set("B = 10000")
		sampleInputText3.set("")
		sampleInputText4.set("")
		sampleInputText5.set("")
	if chooseMethod.get() == 2:
		sampleInputText1.set("n = 108716170")
		sampleInputText2.set("B1 = 1000")
		sampleInputText3.set("B2 = 10000")
		sampleInputText4.set("")
		sampleInputText5.set("")
	if chooseMethod.get() == 3:
		sampleInputText1.set("n = 25147206676")		
		sampleInputText2.set("p_max = 113")
		sampleInputText3.set("B = 10000")
		sampleInputText4.set("B1 = 1000")
		sampleInputText5.set("B2 = 10000")

def nonGUIProcess():
	temp = bytes(linkFile, encoding='utf-8')

	if chooseMethod.get() == 0:
		externalFunctions.m_trival_division(temp)
	elif chooseMethod.get() == 1:
		externalFunctions.m_pollard_rho(temp)
	elif chooseMethod.get() == 2:
		externalFunctions.m_pollard_p1(temp)
	elif chooseMethod.get() == 3:
		externalFunctions.m_full_factorization(temp)

def mainProcess():
	print("Begin process")

	startButton.config(state='disabled')
	loadFileButton.config(state='disabled')
	result.delete('1.0', END)

	t = Thread(target=nonGUIProcess, daemon=True)
	t.start()

	# Hide the main window until the non-GUI stuff is done
	root.withdraw()
	# Create the loading screen
	loading_screen = Toplevel(root)
	loading_screen.title('Loading screen')
	loading_screen.geometry('1200x700')
	loading_label = Label(loading_screen, text="Calculating...")
	loading_label.pack(pady=(150,0))
	pb1 = ttk.Progressbar(loading_screen, orient=HORIZONTAL, length=300, mode='indeterminate')
	pb1.pack(pady=(0,150),expand=True)
	pb1.start(10)

	# While the thread is alive
	while t.is_alive():
		# Update the root so it will keep responding
		root.update()
	
	loading_screen.destroy()
	# Show the main window
	root.deiconify()
	root.focus_force()	

	f = open("result.tmp", "r")
	rs = f.read()
	result.insert(END, rs)
	root.update()
	if os.path.exists("result.tmp"):
  		os.remove("result.tmp")
	
	print("Finish process")
	# startButton.config(state='active')
	loadFileButton.config(state='active')

# ----------------------------Basic parameter-----------------------------
_DIRNAME = os.path.dirname(os.path.abspath(__file__));
so_file = os.path.join(_DIRNAME,"shareLib.so")
externalFunctions = CDLL(so_file)
# -----------------------------------UI-----------------------------------
root = Tk()
root.geometry("1200x700")
root.title("Factorization")
frame = Frame(root, width=1200, height=700)
frame.pack(padx=10)

# A --------Setup area---------
setup_frame = Frame(frame, width=350, height=700)
setup_frame.grid_propagate(False)
setup_frame.grid(row=0, column=0, padx=10, pady=30)

# A-1 Radio group button
methodGroup = LabelFrame(setup_frame, text="Select method", width=300, height=200)
methodGroup.grid(row=0, column=0, sticky="W", padx=5, pady=0, ipadx=0, ipady=0)
chooseMethod = IntVar()
trial_division = Radiobutton(methodGroup, text="Trial division by prime numbers", variable=chooseMethod, value=0, command=changeSampleInput)
trial_division.place(x=10, y=30, anchor="w")
pollard_rho = Radiobutton(methodGroup, text="Pollard's rho", variable=chooseMethod, value=1, command=changeSampleInput)
pollard_rho.place(x=10, y=70, anchor="w")
pollard_p1 = Radiobutton(methodGroup, text="Pollard's p-1", variable=chooseMethod, value=2, command=changeSampleInput)
pollard_p1.place(x=10, y=110, anchor="w")
full_factor = Radiobutton(methodGroup, text="Full factorization", variable=chooseMethod, value=3, command=changeSampleInput)
full_factor.place(x=10, y=150, anchor="w")

# A-2 Load file button
loadFileButton = Button(setup_frame, text="Load input from file", command=select_file, width=34)
loadFileButton.grid(row=1, column=0, sticky="W", padx=5, pady=(20,0), ipadx=0, ipady=0)
inputFilenameText = StringVar()
inputFilenameText.set("Filename: None")
inputFilenameLabel= Label(setup_frame, textvariable= inputFilenameText)
inputFilenameLabel.grid(row=2, column=0, sticky="W", padx=5, pady=0, ipadx=0, ipady=0)
linkFile = ""

# A-3 Loaded input
inputGroup = LabelFrame(setup_frame, text="Your input", width=300, height=180)
inputGroup.grid(row=3, column=0, sticky="W", padx=5, pady=30, ipadx=0, ipady=0)
inputText = []
for i in range(5):
	inputText.append(StringVar())
inputText[0].set('None')
inputLabel1= Label(inputGroup, textvariable= inputText[0])
inputLabel1.place(x=10, y=20, anchor="w")
inputLabel2= Label(inputGroup, textvariable= inputText[1])
inputLabel2.place(x=10, y=50, anchor="w")
inputLabel3= Label(inputGroup, textvariable= inputText[2])
inputLabel3.place(x=10, y=80, anchor="w")
inputLabel4= Label(inputGroup, textvariable= inputText[3])
inputLabel4.place(x=10, y=110, anchor="w")
inputLabel5= Label(inputGroup, textvariable= inputText[4])
inputLabel5.place(x=10, y=140, anchor="w")

# A-4 Start button
startButton = Button(setup_frame, text="Start factorization", command=mainProcess, width=34)
startButton.grid(row=4, column=0, sticky="W", padx=5, pady=40, ipadx=0, ipady=0)
startButton.config(state='disabled')

# A -----------End-------------

# B --------Result area--------
rs_frame = Frame(frame, width=650, height=700)
rs_frame.grid_propagate(False)
rs_frame.grid(row=0, column=1, padx=10, pady=30)

# B-1 Sample input file
sampleGroup = LabelFrame(rs_frame, text="Sample Input File", width=300, height=200)
sampleGroup.grid(row=0, column=0, sticky="W", padx=5, pady=0, ipadx=0, ipady=0)
sampleInputText1 = StringVar()
sampleInputText1.set("n = 2047")
sampleInputLabel1= Label(sampleGroup, textvariable= sampleInputText1)
sampleInputLabel1.place(x=10, y=20, anchor="w")
sampleInputText2 = StringVar()
sampleInputText2.set("p_max = 113")
sampleInputLabel2= Label(sampleGroup, textvariable= sampleInputText2)
sampleInputLabel2.place(x=10, y=50, anchor="w")
sampleInputText3 = StringVar()
sampleInputText3.set("")
sampleInputLabel3= Label(sampleGroup, textvariable= sampleInputText3)
sampleInputLabel3.place(x=10, y=80, anchor="w")
sampleInputText4 = StringVar()
sampleInputText4.set("")
sampleInputLabel4= Label(sampleGroup, textvariable= sampleInputText4)
sampleInputLabel4.place(x=10, y=110, anchor="w")
sampleInputText5 = StringVar()
sampleInputText5.set("")
sampleInputLabel5= Label(sampleGroup, textvariable= sampleInputText5)
sampleInputLabel5.place(x=10, y=140, anchor="w")

# B-2 Result area
rs_label = Label(rs_frame, text="Result").grid(row=1, column=0, padx=5, pady=(20,0))
result = Text(rs_frame, width=75, wrap="word")
result.grid(row=2, column=0, sticky="W", padx=5, pady=0, ipadx=0, ipady=0)
scrollb_y = Scrollbar()
scrollb_y.place(in_=result, relx=1.0, relheight=1.0, bordermode="outside")
scrollb_y.configure(command=result.yview)
# B -----------End-------------

root.mainloop()
# ------------------------------------------------------------------------