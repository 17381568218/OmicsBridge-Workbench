"""Desktop GUI for first-stage metabolomics analysis."""
import os
import queue
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from metabolomics import read_matrix, run_analysis, summarize


DEFAULT_MATRIX = r"G:\softwaropening\cyclicNoteing\output_POS_raw_matrix.csv"


class MetabolomicsApp(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("OmicsBridge Metabolomics")
        self.geometry("980x680")
        self.minsize(860, 560)

        self.log_queue = queue.Queue()
        self.group_values = []
        self.result_dir = ""

        self.matrix_var = tk.StringVar(value=DEFAULT_MATRIX if os.path.exists(DEFAULT_MATRIX) else "")
        self.output_var = tk.StringVar(value=os.path.join(os.getcwd(), "metabolomics_results"))
        self.group_a_var = tk.StringVar()
        self.group_b_var = tk.StringVar()
        self.status_var = tk.StringVar(value="Ready")

        self._build_ui()
        if self.matrix_var.get():
            self.preview_matrix()
        self.after(120, self._drain_log)

    def _build_ui(self):
        root = ttk.Frame(self, padding=12)
        root.pack(fill="both", expand=True)
        root.columnconfigure(0, weight=1)
        root.rowconfigure(2, weight=1)

        header = ttk.Frame(root)
        header.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        ttk.Label(header, text="OmicsBridge Metabolomics", font=("Segoe UI", 17, "bold")).pack(side="left")
        ttk.Label(header, textvariable=self.status_var, foreground="#2563eb").pack(side="right")

        settings = ttk.LabelFrame(root, text="Input and Comparison", padding=10)
        settings.grid(row=1, column=0, sticky="ew", pady=(0, 10))
        settings.columnconfigure(1, weight=1)
        settings.columnconfigure(4, weight=1)

        ttk.Label(settings, text="Matrix CSV").grid(row=0, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(settings, textvariable=self.matrix_var).grid(row=0, column=1, columnspan=4, sticky="ew", pady=4)
        ttk.Button(settings, text="Browse", command=self.browse_matrix).grid(row=0, column=5, padx=(8, 0), pady=4)

        ttk.Label(settings, text="Output").grid(row=1, column=0, sticky="w", padx=(0, 8), pady=4)
        ttk.Entry(settings, textvariable=self.output_var).grid(row=1, column=1, columnspan=4, sticky="ew", pady=4)
        ttk.Button(settings, text="Browse", command=self.browse_output).grid(row=1, column=5, padx=(8, 0), pady=4)

        ttk.Label(settings, text="Group A").grid(row=2, column=0, sticky="w", padx=(0, 8), pady=4)
        self.group_a = ttk.Combobox(settings, textvariable=self.group_a_var, values=self.group_values, width=14)
        self.group_a.grid(row=2, column=1, sticky="w", pady=4)
        ttk.Label(settings, text="Group B").grid(row=2, column=2, sticky="w", padx=(16, 8), pady=4)
        self.group_b = ttk.Combobox(settings, textvariable=self.group_b_var, values=self.group_values, width=14)
        self.group_b.grid(row=2, column=3, sticky="w", pady=4)
        ttk.Button(settings, text="Preview", command=self.preview_matrix).grid(row=2, column=4, sticky="e", padx=(8, 0), pady=4)
        self.run_btn = ttk.Button(settings, text="Run Analysis", command=self.run_clicked)
        self.run_btn.grid(row=2, column=5, sticky="e", padx=(8, 0), pady=4)

        body = ttk.PanedWindow(root, orient="horizontal")
        body.grid(row=2, column=0, sticky="nsew")

        left = ttk.Frame(body, padding=(0, 0, 8, 0))
        right = ttk.Frame(body, padding=(8, 0, 0, 0))
        body.add(left, weight=3)
        body.add(right, weight=2)
        left.rowconfigure(1, weight=1)
        left.columnconfigure(0, weight=1)
        right.rowconfigure(1, weight=1)
        right.columnconfigure(0, weight=1)

        summary_box = ttk.LabelFrame(left, text="Dataset Summary", padding=8)
        summary_box.grid(row=0, column=0, sticky="ew", pady=(0, 10))
        self.summary_text = tk.Text(summary_box, height=8, wrap="word", font=("Consolas", 9))
        self.summary_text.pack(fill="both", expand=True)
        self.summary_text.configure(state="disabled")

        sample_box = ttk.LabelFrame(left, text="Sample Design", padding=8)
        sample_box.grid(row=1, column=0, sticky="nsew")
        sample_box.rowconfigure(0, weight=1)
        sample_box.columnconfigure(0, weight=1)
        self.sample_tree = ttk.Treeview(sample_box, columns=("sample", "mode", "group", "replicate", "qc"), show="headings", height=12)
        for col, width in [("sample", 280), ("mode", 70), ("group", 80), ("replicate", 80), ("qc", 60)]:
            self.sample_tree.heading(col, text=col)
            self.sample_tree.column(col, width=width, anchor="w")
        self.sample_tree.grid(row=0, column=0, sticky="nsew")
        sample_scroll = ttk.Scrollbar(sample_box, orient="vertical", command=self.sample_tree.yview)
        sample_scroll.grid(row=0, column=1, sticky="ns")
        self.sample_tree.configure(yscrollcommand=sample_scroll.set)

        log_box = ttk.LabelFrame(right, text="Run Log", padding=8)
        log_box.grid(row=0, column=0, sticky="nsew")
        right.rowconfigure(0, weight=1)
        self.log_text = tk.Text(log_box, wrap="word", font=("Consolas", 9))
        self.log_text.pack(fill="both", expand=True)
        self.log_text.configure(state="disabled")

        action_box = ttk.LabelFrame(right, text="Results", padding=8)
        action_box.grid(row=1, column=0, sticky="ew", pady=(10, 0))
        ttk.Button(action_box, text="Open Result Folder", command=self.open_results).pack(side="left")
        ttk.Button(action_box, text="Clear Log", command=self.clear_log).pack(side="left", padx=(8, 0))

    def browse_matrix(self):
        path = filedialog.askopenfilename(title="Select metabolomics matrix", filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
        if path:
            self.matrix_var.set(path)
            self.preview_matrix()

    def browse_output(self):
        path = filedialog.askdirectory(title="Select output folder")
        if path:
            self.output_var.set(path)

    def preview_matrix(self):
        path = self.matrix_var.get().strip()
        if not path or not os.path.exists(path):
            messagebox.showwarning("Matrix not found", "Please select a valid matrix CSV.")
            return
        try:
            _, ids, values, _, sample_cols, samples = read_matrix(path)
            summary = summarize(path, ids, values, sample_cols, samples)
            biological_groups = [g["group"] for g in summary["groups"] if not g["group"].upper().startswith("QC")]
            self.group_values = biological_groups
            self.group_a.configure(values=self.group_values)
            self.group_b.configure(values=self.group_values)
            if len(biological_groups) >= 2 and not self.group_a_var.get():
                self.group_a_var.set(biological_groups[0])
                self.group_b_var.set(biological_groups[1])
            self._render_summary(summary)
            self._render_samples(summary["samples"])
            self.status_var.set("Preview loaded")
        except Exception as exc:
            messagebox.showerror("Preview failed", str(exc))

    def _render_summary(self, summary):
        lines = [
            f"Features: {summary['feature_count']}",
            f"Samples: {summary['sample_count']}  Biological: {summary['biological_sample_count']}  QC: {summary['qc_sample_count']}",
            f"Zero: {summary['zero_percent']}%  Missing: {summary['missing_percent']}%",
            "Groups: " + ", ".join(f"{g['group']}({g['count']})" for g in summary["groups"]),
        ]
        if summary["warnings"]:
            lines.append("")
            lines.append("Warnings:")
            lines.extend(f"- {w}" for w in summary["warnings"])
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end")
        self.summary_text.insert("end", "\n".join(lines))
        self.summary_text.configure(state="disabled")

    def _render_samples(self, samples):
        for item in self.sample_tree.get_children():
            self.sample_tree.delete(item)
        for s in samples:
            self.sample_tree.insert("", "end", values=(s["sample_id"], s["mode"], s["group"], s["replicate"], "yes" if s["is_qc"] else "no"))

    def run_clicked(self):
        path = self.matrix_var.get().strip()
        out = self.output_var.get().strip()
        group_a = self.group_a_var.get().strip()
        group_b = self.group_b_var.get().strip()
        if not path or not os.path.exists(path):
            messagebox.showwarning("Matrix not found", "Please select a valid matrix CSV.")
            return
        if not out:
            messagebox.showwarning("Output required", "Please select an output folder.")
            return
        if group_a == group_b:
            messagebox.showwarning("Invalid comparison", "Group A and Group B must be different.")
            return

        self.run_btn.configure(state="disabled")
        self.status_var.set("Running")
        self.clear_log()
        thread = threading.Thread(target=self._run_analysis, args=(path, out, group_a, group_b), daemon=True)
        thread.start()

    def _run_analysis(self, path, out, group_a, group_b):
        try:
            self.log_queue.put("=== Starting metabolomics analysis ===\n")
            self.log_queue.put(f"Input: {path}\n")
            self.log_queue.put(f"Output: {out}\n")
            result = run_analysis(path, out, group_a or None, group_b or None, write_figures=True)
            self.result_dir = out
            self.log_queue.put(f"Features: {result['summary']['feature_count']}\n")
            self.log_queue.put(f"Samples: {result['summary']['sample_count']}\n")
            self.log_queue.put(f"Comparison: {result['comparison']['group_a']} vs {result['comparison']['group_b']}\n")
            self.log_queue.put(f"Significant features: {result['significant_count']}\n")
            self.log_queue.put("Files:\n")
            for filename in result["files"]:
                self.log_queue.put(f"  {filename}\n")
            for warning in result["summary"]["warnings"]:
                self.log_queue.put(f"Warning: {warning}\n")
            self.log_queue.put("=== Done ===\n")
            self.after(0, lambda: self.status_var.set("Done"))
            self.after(0, lambda: messagebox.showinfo("Analysis complete", f"Results saved to:\n{out}"))
        except Exception as exc:
            self.log_queue.put(f"ERROR: {exc}\n")
            self.after(0, lambda: self.status_var.set("Failed"))
            self.after(0, lambda: messagebox.showerror("Analysis failed", str(exc)))
        finally:
            self.after(0, lambda: self.run_btn.configure(state="normal"))

    def _drain_log(self):
        try:
            while True:
                text = self.log_queue.get_nowait()
                self.log_text.configure(state="normal")
                self.log_text.insert("end", text)
                self.log_text.see("end")
                self.log_text.configure(state="disabled")
        except queue.Empty:
            pass
        self.after(120, self._drain_log)

    def clear_log(self):
        self.log_text.configure(state="normal")
        self.log_text.delete("1.0", "end")
        self.log_text.configure(state="disabled")

    def open_results(self):
        path = self.result_dir or self.output_var.get().strip()
        if not path or not os.path.exists(path):
            messagebox.showwarning("No results", "Result folder does not exist yet.")
            return
        os.startfile(path)


if __name__ == "__main__":
    app = MetabolomicsApp()
    app.mainloop()
