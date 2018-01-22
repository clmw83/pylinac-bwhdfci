# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:19:42 2018

@author: clw34
"""

#Helpful resource: https://gist.github.com/dAnjou/2874714

import os
import matplotlib
matplotlib.use('agg')
from flask import Flask, request, redirect, url_for,make_response
from werkzeug import secure_filename
import BWHLeeds
import logging
import tempfile
import traceback
import base64


UPLOAD_FOLDER = os.path.join('.','LeedsTemp')
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)
    
ALLOWED_EXTENSIONS = set(['dcm'])

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

def allowed_file(filename):
    return '.' in filename and \
           filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

@app.route("/", methods=['GET', 'POST'])
def index():
    logtext=""
    if request.method == 'POST':
        
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            tempdir=tempfile.mkdtemp(dir=app.config['UPLOAD_FOLDER'])
            savename=os.path.join(tempdir, filename)
            file.save(savename)
            try:
                leeds=BWHLeeds.BWHLeeds(savename)
                leeds.analyze()
                pdfname = os.path.join(tempdir,"Leeds.pdf")
                leeds.publish_pdf(pdfname)
                data_uri = base64.b64encode(open(pdfname, "rb").read()).decode('ascii')
                logtext+="Successfully processed!\n"
                logtext+="Copy the below lines into the excel workbook:\n"
                logtext+=leeds.excel_lines()
                logtext+="\n"
                logtext+='<button type="button" onclick=DownloadURI("%s")>Download PDF</button>'%data_uri
            except:
                logtext+="Failed processing leeds!"
                logtext+=traceback.format_exc()
            finally:
                os.remove(savename)
                os.remove(pdfname)
                os.rmdir(tempdir)
    return """
    <!doctype html>
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <title>Leeds Phantom Processing</title>
    <script type="text/javascript">
        function DownloadURI(data){
            var blob = new Blob([atob(data)],{type: 'application/pdf'});
            if (window.navigator.msSaveBlob) { // IE
               window.navigator.msSaveOrOpenBlob(blob,"Leeds.pdf");
            }
            else {
                var a = window.document.createElement("a");
                a.href = window.URL.createObjectURL(blob, { type: "application/pdf" });
                a.download = "Leeds.pdf";
                document.body.appendChild(a);
                a.click();
                document.body.removeChild(a);
            }
        }
    </script>
    <h1>Upload new Leeds Phantom File</h1>
    <form action="" method=post enctype=multipart/form-data>
      <p><input type=file name=file>
         <input type=submit value=Process>
    </form>
    <p>%s</p>
    """ % (logtext.replace('\n','<br/>'))

if __name__ == "__main__":
    app.run(host="0.0.0.0",port=5001, debug=True)