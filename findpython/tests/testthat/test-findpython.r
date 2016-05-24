# Copyright (c) 2014 Paul Gilbert
# Copyright (c) 2014 Trevor L. Davis
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

cmdExists <- can_find_python_cmd()

if (cmdExists) {
    context("Any python")
    test_that("python works as expected", {
        cmd <- attr(cmdExists, 'python_cmd')
        # 'sys' should always be installed
        expect_true(can_find_python_cmd(required_modules ='sys'))
        expect_true(is_python_sufficient(cmd, required_modules = 'sys'))

        # 'xxxyyyzzz' is being found and should not exist.")
        expect_false(can_find_python_cmd(required_modules ='xxxyyyzzz', silent=TRUE))
        # Can't figure out how to capture try output via testthat...
        # expect_output(can_find_python_cmd(required_modules ='xxxyyyzzz'), 
        #               "Couldn't find a sufficient Python binary.")
        expect_false(is_python_sufficient(cmd, required_modules ='xxxyyyzzz'))
    })
    python2_exists <- can_find_python_cmd(minimum_version = '2.0', maximum_version = '2.7', silent=TRUE)
    if (python2_exists) {
        cmd <- attr(cmdExists, 'python_cmd')
        context("Python 2")
        test_that("python works as expected", {
            cmd <- attr(python2_exists, 'python_cmd')
            # 'sys' should always be installed
            expect_true(is_python_sufficient(cmd, required_modules = 'sys'))

            # 'xxxyyyzzz' is being found and should not exist.")
            expect_false(is_python_sufficient(cmd, required_modules ='xxxyyyzzz'))
        })
    } else {
        print("Python 2 not found")
    }
    python3_exists <- can_find_python_cmd(minimum_version = '3.0', silent=TRUE)
    if (python3_exists) {
        cmd <- attr(cmdExists, 'python_cmd')
        context("Python 3")
        test_that("python works as expected", {
            cmd <- attr(python3_exists, 'python_cmd')
            # 'sys' should always be installed
            expect_true(is_python_sufficient(cmd, required_modules = 'sys'))

            # 'xxxyyyzzz' is being found and should not exist.")
            expect_false(is_python_sufficient(cmd, required_modules ='xxxyyyzzz'))
        })
    } else {
        print("Python 3 not found")
    }
    # Check Python 4 doesn't exist when testing other machines
    python4_exists <- can_find_python_cmd(minimum_version = '4.0', silent=TRUE)
    if (python4_exists) {
        cmd <- attr(cmdExists, 'python_cmd')
        context("Python 4")
        test_that("python works as expected", {
            cmd <- attr(python4_exists, 'python_cmd')
            # 'sys' should always be installed
            expect_true(is_python_sufficient(cmd, required_modules = 'sys'))

            # 'xxxyyyzzz' is being found and should not exist.")
            expect_false(is_python_sufficient(cmd, required_modules ='xxxyyyzzz'))
        })
    } else {
        print("Python 4 not found")
    }
} else {
    warning('python was not found. No other checks performed.') 
}
