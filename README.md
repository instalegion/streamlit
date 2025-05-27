# Streamlit + Docker scRNA-seq App

Αυτό το repository περιέχει ένα Streamlit app για ανάλυση scRNA-seq δεδομένων με χρήση Scanpy, το οποίο είναι πλήρως Dockerized.

##  Λειτουργίες
- Εισαγωγή αρχείων `.h5ad`
- Προεπεξεργασία δεδομένων (normalization, PCA, clustering)
- DEG (Διαφορική Έκφραση) με πολλαπλές μεθόδους
- Volcano plot, heatmap, dotplot, violin plot
- Αποθήκευση εικόνων
- Docker-ready

##  Εκκίνηση με Docker

```bash
# Βήμα 1: Κατασκευή image
docker build -t my-streamlit-app .

# Βήμα 2: Εκκίνηση container
docker run -p 8501:8501 my-streamlit-app
```

Άνοιξε [http://localhost:8501](http://localhost:8501) για να δεις την εφαρμογή.

##  Requirements

Αν θες να τρέξεις τοπικά:

```bash
pip install -r requirements.txt
streamlit run app.py
```

##  Dependencies
- scanpy
- matplotlib
- pandas
- seaborn
- scanorama
- decoupler
- igraph, leidenalg
- streamlit

##  Συντελεστές
- Σουλιώτης Κωνσταντίνος: UI/UX, Streamlit logic
- Σοφικίτης Αντώνης: DEG ανάλυση, Dockerization

---

*Αναπτύχθηκε για project scRNA-seq με δυνατότητα οπτικοποίησης & επεξεργασίας.*
