import pandas as pd

new_dict = {"DOY": [], "A1": [], "A2": [], "A3": [], "A4": [], "B1": [], "B2": [], "B3": [], "B4": []}
for DOY in range(1, 366):
    with open("GPS navigation messages 2014/brdc" + str(DOY) + "0.14n", "r", encoding = "UTF-8") as input_file:
        head = [next(input_file) for _ in range(5)]
    alphas = head[3].strip()
    while "  " in alphas:
        alphas = alphas.replace("  ", " ")
    betas = head[4].strip()
    while "  " in betas:
        betas = betas.replace("  ", " ")
    alphas = alphas.split(" ")
    betas = betas.split(" ")
    new_dict["DOY"].append(DOY)
    for i in range(4):
        new_dict["A" + str(i + 1)].append(float(alphas[i].replace("D", "E")))
        new_dict["B" + str(i + 1)].append(float(betas[i].replace("D", "E")))
new_df = pd.DataFrame(new_dict)
new_df.to_csv("Klobuchar.csv")