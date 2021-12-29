# -*- coding: utf-8 -*-
####################################################################################################
#
# EntrezDownloader
# Author: Leon Kuchenbecker
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

# The author gave the Frame
# https://github.com/444thLiao thliao write the esummary esearch elink part
# "modified by hychen at 20201019"

import requests
import io
import threading
import time

from concurrent.futures import ThreadPoolExecutor
from concurrent import futures
import pandas as pd
from bs4 import BeautifulSoup


class RequestLimiter:

    def __init__(self, min_wait=0.4):
        """The RequestLimiter class provides functionality to limit the rate at which new requests are made."""
        self.lock = threading.Lock()
        self.last_request = None
        self.min_wait = min_wait

    def wait(self):
        """The wait() function blocks until a minimum wait time from the previous invocation has passed. Thread safe."""
        with self.lock:
            # This is the first request
            if not self.last_request:
                self.last_request = time.time()
                return

            # This is not the first request
            diff = time.time() - self.last_request
            if diff < self.min_wait:
                tsleep = self.min_wait - diff
                time.sleep(tsleep)
            self.last_request = time.time()


class ResultCollector:

    def __init__(self, pbar=None):
        """The ResultCollector class provides functionality for threads to deliver their results."""
        self.pbar = pbar
        self.results = []
        self.failed = []
        self.lock = threading.Lock()

    def add_results(self, results):
        """Adds results to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.results += results
            if self.pbar:
                self.pbar.update(len(results))

    def add_failed(self, ids):
        """Adds failed IDs to the collector. If a progress bar was provided, it updates the progress bar."""
        with self.lock:
            self.failed += ids
            if self.pbar:
                self.pbar.update(len(ids))


class ResultCollectorSearch(ResultCollector):

    def add_results(self, results, ids):
        with self.lock:
            self.results += results
            if self.pbar:
                self.pbar.update(len(ids))


class EntrezDownloader:

    def __init__(self, num_threads=30, batch_size=10, email=None, api_key=None, pbar=False):
        """The EntrezDownloader class enables parallel downloads via the NCBI Entrez interface"""
        self.baseurl = r"https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        self.num_threads = num_threads
        self.batch_size = batch_size
        self.email = email
        self.api_key = api_key
        self.request_limiter = RequestLimiter(
            min_wait=0.4 if not api_key else 0.2)
        self.print_lock = threading.Lock()
        self.pbar = pbar

    def _general_batch(self, db, ids, result_collector, result_func, emode,**kwargs):

        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'db': db,
        }
        if emode == 'esummary':
            if isinstance(ids, str):
                ids = [ids]
            post_data['id'] = ','.join(ids)
        elif emode == 'esearch':
            post_data['term'] = ids

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/{emode}.fcgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    if emode=='esearch' and ' OR ' not in ids:
                        result_collector.add_results(list(zip([_.strip() 
                                                               for _ in ids.split(' OR ')], 
                                                               results)), ids.split(' OR '))
                    else:
                        result_collector.add_results(results,ids.split(' OR '))
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _efetch_batch(self, db, ids, result_collector, result_func, retmode, retype, **kwargs):
        if not retmode:
            retmode = 'text'
        if not retype:
            retype = 'gb'
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ','.join(list(map(str, ids))),
            'db': db,
            'retmode': retmode,
            'rettype': retype
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/efetch.fcgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _esummary_batch(self, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ','.join(ids),
            'db': db,
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/esummary.fcgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def _elink_batch(self, dbfrom, db, ids, result_collector, result_func, **kwargs):
        post_data = {
            'email': self.email,
            'api_key': self.api_key,
            'id': ids,
            'db': db,
            'dbfrom': dbfrom
        }

        post_data.update(kwargs)

        if self.email:
            post_data.update({'email': self.email})

        if self.api_key:
            post_data.update({'api_key': self.api_key})

        error = None
        for i in range(3):  # Retry three times
            try:
                self.request_limiter.wait()
                response = requests.post(
                    f'{self.baseurl}/elink.fcgi', post_data)
                if response.status_code == 200:
                    results = result_func(response.text)
                    result_collector.add_results(results)
                    error = None
                    break
                else:
                    error = f'[STATUS {response.status_code}] An error occurred, you may see a response text here: {response.text}'
            except Exception as e:
                error = f'[UNKNOWN ERROR] {e}'

        if error:
            result_collector.add_failed(ids)
            print(error)

    def elink(self, dbfrom, db, ids,batch_size=None, result_func=lambda x: [x], **kwargs):
        """Interface to the elink database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)
        if not batch_size:
            batch_size = self.batch_size
        fs = []
        for start in range(0, len(ids), self.batch_size):
            num = len(ids)-start
            num = self.batch_size if num > self.batch_size else num
            f = executor.submit(self._elink_batch,
                                db=db,
                                dbfrom=dbfrom,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed

    def efetch(self, db, ids, retmode, retype, result_func=lambda x: [x], batch_size=20, **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        if not batch_size:
            batch_size = self.batch_size
        for start in range(0, len(ids), batch_size):
            num = len(ids)-start
            num = batch_size if num > batch_size else num
            f = executor.submit(self._efetch_batch,
                                db=db,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                retmode=retmode,
                                retype=retype,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed

    def esummary(self, db, ids, batch_size=None,result_func=lambda x: [x], no_pbar= True, **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """

        if self.pbar and not no_pbar:
            from tqdm import tqdm
            results = ResultCollector(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)

        fs = []
        if not batch_size:
            batch_size = self.batch_size
        for start in range(0, len(ids), batch_size):
            num = len(ids)-start
            num = batch_size if num > batch_size else num
            f = executor.submit(self._general_batch,
                                db=db,
                                ids=ids[start:start+num],
                                result_collector=results,
                                result_func=result_func,
                                emode='esummary',
                                **kwargs)
            fs.append(f)

        futures.wait(fs)
        if self.pbar and not no_pbar:
            results.pbar.close()
        return results.results, results.failed

    def esearch(self, db, ids, retmax=0, batch_size=None,result_func=lambda x: [x], **kwargs):
        """Interface to the efetch database.
        result_func: A function to be applied to the response. Must return an iterable.
        """
        if isinstance(ids,str) and ',' in ids:
            ids = [_.strip() for _ in ids.split(',') if _]
        elif isinstance(ids,str):
            ids = [ids]
        if self.pbar:
            from tqdm import tqdm
            results = ResultCollectorSearch(
                pbar=tqdm(total=len(ids), unit='records'))
        else:
            results = ResultCollector()

        executor = ThreadPoolExecutor(max_workers=self.num_threads)
        if not batch_size:
            batch_size = self.batch_size
        if not retmax:
            retmax = self.batch_size*2
        fs = []
        for start in range(0, len(ids), batch_size):
            num = len(ids)-start
            num = batch_size if num > batch_size else num
            f = executor.submit(self._general_batch,
                                db=db,
                                ids=' OR '.join(ids[start:start+num]),
                                result_collector=results,
                                result_func=result_func,
                                emode='esearch',
                                RetMax= retmax,
                                **kwargs)
            fs.append(f)

        futures.wait(fs)
        if self.pbar:
            results.pbar.close()
        return results.results, results.failed


if __name__ == "__main__":
    edl = EntrezDownloader(
        # An email address. You might get blocked by the NCBI without specifying one.
        email="eacochen@163.com",
        # An API key. You can obtain one by creating an NCBI account. Speeds things up.
        api_key='df1607793b0171d9f7fd39585f90a3491007',
        num_threads=30,                       # The number of parallel requests to make
        batch_size=500,                        # The number of IDs to fetch per request
        pbar=True                             # Enables a progress bar, requires tqdm package
    )
