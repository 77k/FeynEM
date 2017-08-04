/***
 * simple httpd for REST style API
 * DON'T EXPOSE TO THE INTERWEBz
 * snaske@77k.eu
 */
#ifndef NAIVE_HTTPD_H
#define NAIVE_HTTPD_H
#include <string>
#include <vector>
#include <limits>

namespace naive_httpd {

    enum limits
    {
        MAX_METHOD_LENGTH = 8,
        MAX_FT_LENGTH = 8 
    };

    enum status_type
    {
        ok = 200,
        created = 201,
        accepted = 202,
        no_content = 204,
        multiple_choices = 300,
        moved_permanently = 301,
        moved_temporarily = 302,
        not_modified = 304,
        bad_request = 400,
        unauthorized = 401,
        forbidden = 403,
        not_found = 404,
        internal_server_error = 500,
        not_implemented = 501,
        bad_gateway = 502,
        service_unavailable = 503
    };


    struct ext2type
    {
        std::map< std::string, std::string > _f2ms;
        ext2type(){
            _f2ms.insert(std::make_pair("gif", "image/gif"));
            _f2ms.insert(std::make_pair("htm", "text/html"));
            _f2ms.insert(std::make_pair("html", "text/html"));
            _f2ms.insert(std::make_pair("css", "text/css"));
            _f2ms.insert(std::make_pair("jpg", "image/jpeg"));
            _f2ms.insert(std::make_pair("png", "image/png"));
            _f2ms.insert(std::make_pair("json", "application/json"));
            _f2ms.insert(std::make_pair("svg", "image/svg+xml"));
        };
        std::string operator [] (const std::string fe)
        {

            if(([&](){ return (fe.length() > limits::MAX_FT_LENGTH) ?
                        _f2ms.end() : _f2ms.find(fe);})() == _f2ms.end())
                return "plain/text";
            return _f2ms[fe];
        }
    };

    struct header
    {
        std::string name;
        std::string value;
    };

    struct request
    {
        std::string method;
        std::string uri;
        float http_version;
        int http_version_major;
        int http_version_minor;
        std::vector<header> headers;
        std::string content;
    };

    const char http_id[] = "HTTP/";

    class req_parser : boost::asio::coroutine
    {
        public:
#include <boost/asio/yield.hpp>
            template < typename It >
		    boost::tuple < boost::tribool, It > parse(request &req,
				    It beg, It end)
                {
                    It beg0, end0;
                    std::string request_string;
                    beg0 = beg; end0 = end;
                    boost::tribool res;

                    std::string le_string = std::string(beg0, end0);
                    std::string space0_str;
                    std::string space1_str;
                    std::string http_str;
                    std::string crlf0_str;
                    boost::regex meth_exp {"^([A-Z]{1,8})?( )?(.*)"};
                    boost::smatch matches;


                    boost::regex_match(le_string, matches, meth_exp);
                    req.method.clear();
                    req.method = matches[1];
                    if(req.method.empty()) return boost::make_tuple(false, beg);
		    space0_str = matches[2];
		    if(space0_str.empty()) return boost::make_tuple(
				    static_cast<boost::tribool>(boost::indeterminate), beg);
		    beg += req.method.size() + space0_str.size();
		    le_string = matches[3];
		    boost::regex http_exp{
			    "^([/A-Za-z0-9_\\+%\\.-]+)?( )?(HTTP/(\\d\\.\\d))?(\\r\\n)?(.*)"};
		    boost::regex_match(le_string, matches, http_exp);
		    req.uri = matches[1];
                    beg += req.uri.size();
                    if(req.uri.empty()) return boost::make_tuple(false, beg);
                    space1_str = matches[2];
                    beg += space1_str.size();
                    if(space1_str.empty()) return boost::make_tuple(false, beg);
                    http_str = matches[3];
                    beg += http_str.size();
                    if(http_str.empty()) return boost::make_tuple(false, beg);
                    crlf0_str = matches[5];
                    beg += crlf0_str.size();
                    if(crlf0_str.empty()) return boost::make_tuple(false, beg);

                    while(beg != end)
                    {
                        res = process_headers(req, *beg++);
                        if(res || !res) return boost::make_tuple(res, beg);
                    }

                    res = boost::indeterminate;
                    return boost::make_tuple(res, beg);
                }
    
#include <boost/asio/unyield.hpp>

    std::string req_str;
    std::size_t content_length_;
    boost::tribool process_headers(request& req, char input);
    static bool is_char(int c); static bool is_tspecial(int c); //replace that foo
};

struct response
{
    /// The status of the response.
    status_type status;
    std::string bla;
    std::shared_ptr< std::map< status_type, std::string > > status2string_p;

    response()
    {
        if(!status2string_p)
        {
            status2string_p.reset( new std::map< status_type, std::string >());
            status2string_p->insert(std::make_pair(status_type::ok, "HTTP/1.1 200 OK\r\n"));
            status2string_p->insert(std::make_pair(status_type::created, "HTTP/1.1 201 Created\r\n"));
            status2string_p->insert(std::make_pair(status_type::accepted, "HTTP/1.1 202 Accepted\r\n"));
            status2string_p->insert(std::make_pair(status_type::no_content, "HTTP/1.1 204 No Content\r\n"));
            status2string_p->insert(std::make_pair(status_type::multiple_choices, "HTTP/1.1 300 Muliple Coices\r\n"));
            status2string_p->insert(std::make_pair(status_type::moved_permanently, "HTTP/1.1 301 Moved Permanently\r\n"));
            status2string_p->insert(std::make_pair(status_type::moved_temporarily, "HTTP/1.1 302 Moved Temporarily\r\n"));
            status2string_p->insert(std::make_pair(status_type::not_modified, "HTTP/1.1 304 Not Modified\r\n"));
            status2string_p->insert(std::make_pair(status_type::bad_request, "HTTP/1.1 400 Bad Request\r\n"));
            status2string_p->insert(std::make_pair(status_type::unauthorized, "HTTP/1.1 401 Unauthorized\r\n"));
            status2string_p->insert(std::make_pair(status_type::forbidden, "HTTP/1.1 403 Forbidden\r\n"));
            status2string_p->insert(std::make_pair(status_type::not_found, "HTTP/1.1 404 Not Found\r\n"));
            status2string_p->insert(std::make_pair(status_type::internal_server_error, "HTTP/1.1 500 Internal Server Error\r\n"));
            status2string_p->insert(std::make_pair(status_type::not_implemented, "HTTP/1.1 501 Not Implemented\r\n"));
            status2string_p->insert(std::make_pair(status_type::bad_gateway, "HTTP/1.1 502 Bad Gateway\r\n"));
            status2string_p->insert(std::make_pair(status_type::service_unavailable, "HTTP/1.1 503 Service Unavailable\r\n"));
        }

    }

    std::string& operator [] (const status_type s)
    {
        return (*status2string_p)[s];
    }

    std::vector<header> headers;

    std::string content;

    std::vector< boost::asio::const_buffer > to_buffers();
    static response default_response(status_type status);
};

static                    std::string doc_root_; //clean it up fucker

namespace method_strings //what the fuck have i done?
{
    const std::string GET = "GET";
    const std::string HEAD = "HEAD";
    const std::string POST = "POST";
    const std::string PUT = "PUT";
    const std::string DELETE = "DELETE";
    const std::string TRACE = "TRACE";
    const std::string CONNECT = "CONNECT";
}

struct response;

template < class _T >
class REST_handler
{
    public:
        REST_handler(const std::string& doc_root)
        {
            _T::init();
            doc_root_ = doc_root;
            methods_.insert(make_pair(method_strings::GET, GET()));
            methods_.insert(make_pair(method_strings::HEAD, HEAD()));
            methods_.insert(make_pair(method_strings::POST, POST()));
            methods_.insert(make_pair(method_strings::PUT, PUT()));
            methods_.insert(make_pair(method_strings::DELETE, DELETE()));
            methods_.insert(make_pair(method_strings::TRACE, TRACE()));
            methods_.insert(make_pair(method_strings::CONNECT, CONNECT()));
        }

        response& operator()(const request& req, response& resp)
        {

	    //std::cerr << "base req11 = " << req.content << std::endl;
            if(([&](){ return (req.method.length() > limits::MAX_METHOD_LENGTH) ?
                        methods_.end() : methods_.find(req.method);})() == methods_.end())
            {   
                return resp = response::default_response(bad_request);
            }
	    //std::cerr << "base req = " << req.content << std::endl;
            boost::regex val_0 {"^(([^:/?#]+):)?(//([^/?#]*))?([^?#]*)(\\?([^#]*))?(#(.*))?"};
            boost::smatch matches;
            if(!boost::regex_match(req.uri, matches, val_0))
            {
                return resp = response::default_response(bad_request);
            }
            std::string path(matches[5]);
            if (path.empty() || path[0] != '/' || path.find("..") != std::string::npos)
            {
                return resp = response::default_response(bad_request);
            }

            methods_[req.method](req, resp);
            return resp;
        }

        typedef typename _T::GET GET;
        typedef typename _T::HEAD HEAD;
        typedef typename _T::POST POST;
        typedef typename _T::PUT PUT;
        typedef typename _T::DELETE DELETE;
        typedef typename _T::TRACE TRACE;
        typedef typename _T::CONNECT CONNECT;

    private:
        std::map< decltype(request::method), std::function< void (const request&, response&) > > methods_;
        std::stack< std::string > stack_; 

};


class blarghd : boost::asio::coroutine
{
    public:
        explicit blarghd(boost::asio::io_service& io_service, const std::string& ip_addresss, const std::string& port, boost::function<void(const request&, response&)> req_handler);
        void operator()(
                boost::system::error_code ec = boost::system::error_code(),
                std::size_t length = 0);
    private:
        typedef boost::asio::ip::tcp tcp;
        boost::function<void(const request&, response&)> m_req_handler;
        boost::shared_ptr<tcp::acceptor> m_acceptor;
        boost::shared_ptr<tcp::socket> m_socket;
        boost::shared_ptr<boost::array<char, 8192> > m_buffer;

        boost::shared_ptr<request> m_req;
        boost::shared_ptr<response> m_response;
        req_parser m_req_parser;
        boost::tribool m_valid_req;

};

#include <boost/asio/yield.hpp>

    //first match HTTP... or extract methods?.
    boost::tribool req_parser::process_headers(request& req, char c)
    {
        reenter (this)
        {
            req.http_version_major = 0;
            req.http_version_minor = 0;
            req.headers.clear();
            req.content.clear();
            content_length_ = 0;
            while ((is_char(c) && !iscntrl(c) && !is_tspecial(c) && c != '\r')
                    || (c == ' ' || c == '\t'))
            {
                if (c == ' ' || c == '\t')
                {
                    if (req.headers.empty()) return false;
                    while (c == ' ' || c == '\t')
                        yield return boost::indeterminate;
                }
                else
                {
                    req.headers.push_back(header());

                    while (is_char(c) && !iscntrl(c) && !is_tspecial(c) && c != ':')
                    {
                        req.headers.back().name.push_back(c);
                        yield return boost::indeterminate;
                    }

                    if (c != ':') return false;
                    yield return boost::indeterminate;
                    if (c != ' ') return false;
                    yield return boost::indeterminate;
                }

                while (is_char(c) && !iscntrl(c) && c != '\r')
                {
                    req.headers.back().value.push_back(c);
                    yield return boost::indeterminate;
                }

                if (c != '\r') return false;
                yield return boost::indeterminate;
                if (c != '\n') return false;
                yield return boost::indeterminate;
            }

            if (c != '\r') return false;
            yield return boost::indeterminate;
            if (c != '\n') return false;

            for (std::size_t i = 0; i < req.headers.size(); ++i)
            {
                if(req.headers[i].name == "Content-Length")
                {
                    try
                    {
                        content_length_ = std::stoi(req.headers[i].value);
                    }
                    catch (...)//invalid cast ... out of range
                    {
                        return false;
                    }
                }
            }

            while (req.content.size() < content_length_)
            {
                yield return boost::indeterminate;
                req.content.push_back(c);
            }
        }

        return true;
    }

#include <boost/asio/unyield.hpp>

    bool req_parser::is_char(int c)
    {
        return c >= 0 && c <= 127;
    }


    bool req_parser::is_tspecial(int c)
    {
        switch (c)
        {
            case '(': case ')': case '<': case '>': case '@':
            case ',': case ';': case ':': case '\\': case '"':
            case '/': case '[': case ']': case '?': case '=':
            case '{': case '}': case ' ': case '\t':
                return true;
            default:
                return false;
        }
    }


    const std::string cat = "<pre> yo</pre>";        
    char kvsep[] = { ':', ' ' };
    char crlf[] = { '\r', '\n' };

    response response::default_response(status_type status)
    {

        std::string resp_str;
        if(status != status_type::ok)
        {
            resp_str = "<html><head><title>" + std::to_string(status) + "</title>" + "</head>" + "<body>" + cat + "\n <img src=\"./blargh.gif\"/>" + " </body>" + "</html>";
        }
        response resp;
        resp.status = status;
        resp.content = resp_str;
        resp.headers.resize(2);
        resp.headers[0].name = "Content-Length";
        resp.headers[0].value = boost::lexical_cast<std::string>(resp.content.size());
        resp.headers[1].name = "Content-Type";
        resp.headers[1].value = "text/html";
        return resp;
    }

    std::vector<boost::asio::const_buffer> response::to_buffers()
    {

        std::vector<boost::asio::const_buffer> v_b;
        v_b.push_back(boost::asio::buffer((*this)[status_type::ok]));
        std::for_each(headers.begin(), headers.end(), [&](header &hd)
                { 
                v_b.push_back(boost::asio::buffer(hd.name));
                v_b.push_back(boost::asio::buffer(kvsep));
                v_b.push_back(boost::asio::buffer(hd.value));
                v_b.push_back(boost::asio::buffer(crlf));

                });
        v_b.push_back(boost::asio::buffer(crlf));
        v_b.push_back(boost::asio::buffer(content));
        return v_b;
    }



    blarghd::blarghd(boost::asio::io_service& io_service, const std::string& ip_address, const std::string& port, boost::function<void(const request&, response&)> req_handler)
        : m_req_handler(req_handler)
    {
        tcp::resolver resolver(io_service);
        tcp::resolver::query query(ip_address, port);
        m_acceptor.reset(new tcp::acceptor(io_service, *resolver.resolve(query)));
    }

#include <boost/asio/yield.hpp>
    //for details look into the boost coroutine documentation and consult their examples, we simply use the example code for rapid prototyping
    void blarghd::operator()(boost::system::error_code ec, std::size_t length)
    {
        if (!ec)
        {
            reenter (this)
            {
                do
                {
                    m_socket.reset(new tcp::socket(m_acceptor->get_io_service()));
                    yield m_acceptor->async_accept(*m_socket, *this);
                    fork blarghd(*this)();

                } while (is_parent());

                m_buffer.reset(new boost::array<char, 8192>);
                m_req.reset(new request);

                do
                {
                    yield m_socket->async_read_some(boost::asio::buffer(*m_buffer), *this);

                    boost::tie(m_valid_req, boost::tuples::ignore)
                        = m_req_parser.parse(*m_req,
                                m_buffer->data(), m_buffer->data() + length);

                } while (boost::indeterminate(m_valid_req));

                m_response.reset(new response);

                if (m_valid_req)
                {
                    m_req_handler(*m_req, *m_response);
                }
                else
                {
                    *m_response = response::default_response(status_type::bad_request);
                }

                yield boost::asio::async_write(*m_socket, m_response->to_buffers(), *this);

                m_socket->shutdown(tcp::socket::shutdown_both, ec);
            }
        }

    }

#include <boost/asio/unyield.hpp>




}

#endif

